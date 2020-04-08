# Copyright (C) 1999-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# JNET -z format parsing consists of:
#   ALIGNMENT
#
###########################################################################
package Bio::Parse::Format::JNETZ;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Parse::Record);

#delimit full JNETZ entry
my $JNETZ_START          = '^START PRED';
my $JNETZ_END            = '^END PRED';
my $JNETZ_Null           = '^\s*$';#'

my $JNETZ_ALIGNMENT      = $JNETZ_START;
my $JNETZ_ALIGNMENTend   = $JNETZ_END;


#Consume one entry-worth of input on text stream associated with $file and
#return a new JNETZ instance.
sub get_entry {
    my $text = shift;
    my $line = '';
    my $data = 0;

    while ($text->getline(\$line)) {

        #start of entry
        if ($line =~ /$JNETZ_START/o and !$data) {
            $text->start_count();
            $data = 1;
            next;
        }

        #end of entry
        if ($line =~ /$JNETZ_END/o and $data) {
            $text->stop_count_at_end();
            last;
        }
    }
    return 0  unless $data;

    new Bio::Parse::Format::JNETZ(undef, $text, $text->get_start(), $text->get_stop());
}

#Parse one entry
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #ALIGNMENT lines
        if ($line =~ /$JNETZ_ALIGNMENT/o) {
            $scan->scan_until($JNETZ_ALIGNMENTend);
            $self->push_record('ALIGNMENT',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #blank line or empty record: ignore
        next  if $line =~ /$JNETZ_Null/o;

        #terminal line: ignore
        next  if $line =~ /$JNETZ_END/o;

        #default
        $self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::JNETZ::ALIGNMENT;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #initialise these:
    $self->{'query'} = [];
    $self->{'final'} = [];
    $self->{'conf'}  = [];
    $self->{'align'} = [];

    while (defined ($line = $scan->read_line)) {

        if ($line =~
            /^\s*(\S)
            \s+(\S)
            \s+\|
            \s+(\S)
            \s+.
            \s+(\d)
            /xo
           ) {
            push @{ $self->{'query'}}, $1;
            push @{ $self->{'final'}}, $2;
            push @{ $self->{'align'}}, $3;
            push @{ $self->{'conf'}},  $4;
            next;
        }

        #header line: ignore
        next  if $line =~ /$JNETZ_ALIGNMENT/;

        #blank line or empty record: ignore
        next  if $line =~ /$JNETZ_Null/o;

        #default
        $self->warn("unknown field: $line");
    }

    $self;#->examine;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> %s\n",   'query', scalar $self->get_query;
    $s .= sprintf "$x%20s -> %s\n",   'final', scalar $self->get_final;
    $s .= sprintf "$x%20s -> %s\n",   'conf',  scalar $self->get_conf;
    $s .= sprintf "$x%20s -> %s\n",   'align', scalar $self->get_align;
    return $s;
}

sub get_query {
    my $self = shift;
    my @a = @{$self->{'query'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_final {
    my $self = shift;
    my @a = @{$self->{'final'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_conf {
    my $self = shift;
    my @a = @{$self->{'conf'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_align {
    my $self = shift;
    my @a = @{$self->{'align'}};
    return \@a    if wantarray;
    return join("", @a);
}


###########################################################################
1;
