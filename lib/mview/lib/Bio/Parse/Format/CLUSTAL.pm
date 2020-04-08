# Copyright (C) 1998-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::CLUSTAL;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Parse::Record);


#delimit full CLUSTAL entry
my $CLUSTAL_START          = '^\s*CLUSTAL';
my $CLUSTAL_END            = $CLUSTAL_START;

#CLUSTAL record types
my $CLUSTAL_ALIGNMENT      = '^\s*\S+\s+\S+(?:\s+\d+)?\s*$';
my $CLUSTAL_ALIGNMENTend   = $CLUSTAL_START;
my $CLUSTAL_HEADER         = $CLUSTAL_START;
my $CLUSTAL_HEADERend      = "(?:$CLUSTAL_ALIGNMENT|$CLUSTAL_HEADER)";
my $CLUSTAL_Null           = '^\s*$';#'


#Consume one entry-worth of input on text stream associated with $file and
#return a new CLUSTAL instance.
sub get_entry {
    my $text = shift;
    my $line = '';
    my $data = 0;

    while ($text->getline(\$line)) {

        #start of entry
        if ($line =~ /$CLUSTAL_START/o and !$data) {
            $data = 1;
            next;
        }

        #consume rest of stream
        if ($line =~ /$CLUSTAL_END/o and $data) {
            $text->stop_count_at_start();
            last;
        }
    }
    return 0  unless $data;

    new Bio::Parse::Format::CLUSTAL(undef, $text, $text->get_start(), $text->get_stop());
}

#Parse one entry
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #HEADER lines
        if ($line =~ /$CLUSTAL_HEADER/o) {
            $scan->scan_until($CLUSTAL_HEADERend);
            $self->push_record('HEADER',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #consume data

        #ALIGNMENT lines
        if ($line =~ /$CLUSTAL_ALIGNMENT/o) {
            $scan->scan_until($CLUSTAL_ALIGNMENTend);
            $self->push_record('ALIGNMENT',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$CLUSTAL_Null/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::CLUSTAL::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'version'} = '';
    $self->{'major'}   = '';
    $self->{'minor'}   = '';

    #consume Name lines
    while (defined ($line = $scan->read_line)) {

        #first part of CLUSTAL line
        if ($line =~ /^
            \s*
            CLUSTAL
            \s+
            (([^\(\s]+)    #major version, eg., W
            \s*
            \((\S+)\))     #minor version, eg., 1.70
            /xo) {

            $self->test_args(\$line, $1, $2, $3);
            (
             $self->{'version'},
             $self->{'major'},
             $self->{'minor'},
            ) = ($1, $2, $3);

        }

        #first part of T-COFFEE line
        if ($line =~ /^
            \s*
            CLUSTAL\sFORMAT\sfor\s((T-COFFEE)
            \s+
            Version_(\S+)),     #version number, eg., 1.37
            /xo) {

            $self->test_args(\$line, $1, $2, $3);
            (
             $self->{'version'},
             $self->{'major'},
             $self->{'minor'},
            ) = ($1, $2, $3);

        }

        #ignore any other text
    }

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> %s\n",   'version', $self->{'version'};
    $s .= sprintf "$x%20s -> %s\n",   'major',   $self->{'major'};
    $s .= sprintf "$x%20s -> %s\n",   'minor',   $self->{'minor'};
    return $s;
}


###########################################################################
package Bio::Parse::Format::CLUSTAL::ALIGNMENT;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'id'}    = [];
    $self->{'seq'}   = {};
    $self->{'match'} = '';

    my $off = 0;

    while (defined ($line = $scan->read_line(1))) {

        #match symbols, but only if expected
        if ($off and $line !~ /[^*:. ]/) {
            local $^W=0;
            $line = substr($line, $off);
            $self->{'match'} .= $line;
            $off = 0;
            next;
        }

        #id/sequence (/optional number)
        if ($line =~ /^\s*(\S+)\s+(\S+)(?:\s+\d+)?$/o) {
            $self->test_args(\$line, $1, $2);
            push @{$self->{'id'}}, $1    unless exists $self->{'seq'}->{$1};
            $self->{'seq'}->{$1} .= $2;
            $off = length($line) - length($2);
            next;
        }

        next    if $line =~ /$CLUSTAL_Null/o;

        #default
        $self->warn("unknown field: $line");
    }

    #line length check (ignore 'match' as this may be missing)
    if (defined $self->{'id'}->[0]) {
        $off = length $self->{'seq'}->{$self->{'id'}->[0]};
        foreach my $id (keys %{$self->{'seq'}}) {
            $line = $self->{'seq'}->{$id};
            my $len = length $line;
            #warn "$off, $len, $id\n";
            if ($len != $off) {
                $self->die("length mismatch for '$id' (expect $off, saw $len):\noffending sequence: [$line]\n");
            }
        }
    }

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $i (@{$self->{'id'}}) {
        $s .= sprintf "$x%20s -> %-15s %s\n", 'seq', $i, $self->{'seq'}->{$i};
    }
    $s .= sprintf "$x%20s -> %-15s %s\n", 'match', '', $self->{'match'};
    return $s;
}


###########################################################################
1;
