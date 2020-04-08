# Copyright (C) 1998-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::MSF;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Parse::Record);


#delimit full MSF entry
#my $MSF_START          = '^\s*(?:\S+)?\s+MSF:';
#my $MSF_START          = '^(?:PileUp|\s*(?:\S+)?\s+MSF:)';
my $MSF_START          = '^(?:PileUp|.*\s+MSF:)';
my $MSF_END            = '^PileUp';

#MSF record types
my $MSF_HEADER         = $MSF_START;
my $MSF_HEADERend      = '^\s*Name:';
my $MSF_NAME           = $MSF_HEADERend;
my $MSF_NAMEend        = '^\/\/';
my $MSF_ALIGNMENT      = '^\s*\S+\s+\S';
my $MSF_ALIGNMENTend   = $MSF_START;
my $MSF_Null           = '^\s*$';#'


#Consume one entry-worth of input on text stream associated with $file and
#return a new MSF instance.
sub get_entry {
    my $text = shift;
    my $line = '';
    my $data = 0;

    while ($text->getline(\$line)) {

        #start of entry
        if ($line =~ /$MSF_START/o and !$data) {
            $text->start_count();
            $data = 1;
            next;
        }

        #consume rest of stream
        if ($line =~ /$MSF_END/o and $data) {
            $text->stop_count_at_start();
            last;
        }
    }
    return 0  unless $data;

    new Bio::Parse::Format::MSF(undef, $text, $text->get_start(), $text->get_stop());
}

#Parse one entry
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #HEADER lines
        if ($line =~ /$MSF_HEADER/o) {
            $scan->scan_until($MSF_HEADERend);
            $self->push_record('HEADER',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #consume data

        #NAME lines
        if ($line =~ /$MSF_NAME/o) {
            $scan->scan_until($MSF_NAMEend);
            $self->push_record('NAME',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #ALIGNMENT lines
        if ($line =~ /$MSF_ALIGNMENT/o) {
            $scan->scan_until($MSF_ALIGNMENTend);
            $self->push_record('ALIGNMENT',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$MSF_Null/o;

        #end of NAME section: ignore
        next    if $line =~ /$MSF_NAMEend/o;

        #default
        $self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::MSF::HEADER;

use Bio::Parse::Strings qw(strip_trailing_space);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #consume Name lines
    while (defined ($line = $scan->read_line)) {

        #MSF line
        if ($line =~ /^
            \s*
            ((?:.+)?)
            MSF\:\s+(\d+)
            \s+
            Type\:\s+(\S+)
            \s*
            ((?:.+)?)
            Check\:\s+(\d+)
            \s+\.\.
            /xo) {

            $self->test_args(\$line, $2,$3,$5);
            (
             $self->{'file'},
             $self->{'msf'},
             $self->{'type'},
             $self->{'data'},
             $self->{'check'},
            ) = (strip_trailing_space($1),
                 $2,$3,
                 strip_trailing_space($4),
                 $5);
        }

        #ignore any other text
    }

    $self->warn("missing MSF data\n")  unless exists $self->{'msf'};

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> '%s'\n", 'file',   $self->{'file'};
    $s .= sprintf "$x%20s -> %s\n",   'msf',    $self->{'msf'};
    $s .= sprintf "$x%20s -> %s\n",   'type',   $self->{'type'};
    $s .= sprintf "$x%20s -> '%s'\n", 'data',   $self->{'data'};
    $s .= sprintf "$x%20s -> %s\n",   'check',  $self->{'check'};
    return $s;
}


###########################################################################
package Bio::Parse::Format::MSF::NAME;

use Bio::Parse::Strings qw(strip_trailing_space);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'seq'}   = {};
    $self->{'order'} = [];

    #consume Name lines
    while (defined ($line = $scan->read_line)) {
        my $id = "";

        if ($line =~ /^(\s*Name:\s+(.+))Len:\s+\d/) {
            $line = substr($line, length($1));
            $id = $2;
            if ($id =~ /^(.*)\s+oo\s*$/) { #weird clustal insertion
                $id = $1;
            }
            $id = strip_trailing_space($id);
            #warn "[$id] [$line]\n";

            if ($line =~ /^
                Len\:\s+(\d+)     #sequence length
                \s+
                Check\:\s+(\S+)   #checksum
                \s+
                Weight\:\s+(\S+)  #sequence weight
                /xo) {

                $self->test_args(\$line, $1,$2,$3);
                $self->{'seq'}->{$id} = {
                    'length' => $1,
                    'check'  => $2,
                    'weight' => $3,
                };
                push @{$self->{'order'}}, $id;
                next;
            }
        }

        next  if $line =~ /$MSF_Null/;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $i (@{$self->{'order'}}) {
        $s .= sprintf "$x%20s -> %-15s %s=%5s %s=%5s %s=%5s\n",
            'seq',    $i,
            'length', $self->{'seq'}->{$i}->{'length'},
            'check',  $self->{'seq'}->{$i}->{'check'},
            'weight', $self->{'seq'}->{$i}->{'weight'};
    }
    return $s;
}


###########################################################################
package Bio::Parse::Format::MSF::ALIGNMENT;

use Bio::Parse::Strings qw(strip_leading_space strip_trailing_space);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'seq'} = {};

    my $maxnamelen = 0;
    my $sibling = $self->get_parent(1)->{'blockkeeper'}->get_block('NAME')->{'record'};

    foreach my $id ( @{$sibling->{'order'}} ) {
        my $len = length($id);
        $maxnamelen = $len  if $len > $maxnamelen;
    }
    #warn $maxnamelen;

    while (defined ($line = $scan->read_line)) {

        #start/end positions
        next  if $line =~ /^\s*\d+\s+\d+$/o;

        #end position
        next  if $line =~ /^\s*\d+\s*$/o;

        #id/sequence
        if ($line =~ /^\s*(.{$maxnamelen})\s+(.*)$/o) {
            my $id;
            $id = strip_leading_space($1);
            $id = strip_trailing_space($id);
            $self->test_args(\$line, $id, $2);
            $self->{'seq'}->{$id} .= $2;
            next;
        }

        next  if $line =~ /$MSF_Null/;

        #default
        $self->warn("unknown field: $line");
    }

    foreach (keys %{$self->{'seq'}}) {
        $self->{'seq'}->{$_} =~ s/ //g;
    }

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $i (sort keys %{$self->{'seq'}}) {
        $s .= sprintf "$x%20s -> %-15s =  %s\n", 'seq', $i, $self->{'seq'}->{$i};
    }
    return $s;
}


###########################################################################
1;
