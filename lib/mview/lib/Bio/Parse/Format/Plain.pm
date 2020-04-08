# Copyright (C) 1998-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
#
# The simplest input alignment is a column of id's and a column of aligned
# sequences of identical length with the entire alignment in one block.
#
###########################################################################
package Bio::Parse::Format::Plain;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

my $BLANK           = '^\s*$';
my $COMMENT         = '^\s*\#';

#subrecords
my $ALIGNMENT_BEGIN = '^\s*\S+\s+\S';
my $ALIGNMENT_GROUP = '^\s*(\S+)\s+(\S+)\s*$';
my $ALIGNMENT_END   = $BLANK;

#full record
my $DATA_BEGIN      = $ALIGNMENT_BEGIN;
my $DATA_END        = $ALIGNMENT_END;

my $DEBUG = 0;

#Consume one entry-worth of input on text stream associated with $file and
#return a new Plain instance.
sub get_entry {
    my $text = shift;
    my $line = '';
    my $data = 0;

    #warn "BEGIN ENTRY\n"  if $DEBUG;

    while ($text->getline(\$line)) {

        #initial blank
        if ($line =~ /$BLANK/o and !$data) {
            #warn " LEADING BLANK\n"  if $DEBUG;
            next;
        }

        #initial comment
        if ($line =~ /$COMMENT/o and !$data) {
            #warn " LEADING COMMENT\n"  if $DEBUG;
            next;
        }

        #start of data
        if ($line =~ /$DATA_BEGIN/o and !$data) {
            #warn " START OF DATA\n"  if $DEBUG;
            $text->start_count();
            $data = 1;
            next;
        }

        #end of data
        if ($line =~ /$DATA_END/o and $data) {
            #warn " END OF DATA\n"  if $DEBUG;
            $text->stop_count_at_start();
            last;
        }
    }

    #warn "END ENTRY\n"  if $DEBUG;

    return 0  unless $data;

    new Bio::Parse::Format::Plain(undef, $text, $text->get_start(), $text->get_stop());
}

#Parse one entry
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #warn "BEGIN PARSE\n"  if $DEBUG;

    while (defined ($line = $scan->read_line)) {

        #alignment line; can include sequence identifiers prefixed with '#'
        #for special colouring treatment
        if ($line =~ /$ALIGNMENT_BEGIN/o) {

            $scan->scan_until($ALIGNMENT_END);
            $self->push_record('ALIGNMENT',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #comment: not-matching 'id seq' pattern above; test last
        next  if $line =~ /$COMMENT/o;

        #blank: ignore; test last
        next  if $line =~ /$BLANK/o;

        $self->warn("unknown field: $line");
    }

    #warn "END PARSE\n"  if $DEBUG;

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::Plain::ALIGNMENT;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'id'}    = [];
    $self->{'seq'}   = {};

    #warn "BEGIN PARSE ALIGNMENT\n"  if $DEBUG;

    while (defined ($line = $scan->read_line(1))) {

        #id sequence
        if ($line =~ /$ALIGNMENT_GROUP/o) {
            #warn "ALIGNMENT_GROUP\n"  if $DEBUG;
            $self->test_args(\$line, $1, $2);
            push @{$self->{'id'}}, $1    unless exists $self->{'seq'}->{$1};
            $self->{'seq'}->{$1} .= $2;
            next;
        }

        #comment: not-matching 'id seq' pattern above; test last
        next  if $line =~ /$COMMENT/o;

        #default
        $self->warn("unknown field: $line");
    }

    #warn "END PARSE ALIGNMENT\n"  if $DEBUG;

    #line length check
    my $off = 0;
    if (defined $self->{'id'}->[0]) {
        $off = length $self->{'seq'}->{$self->{'id'}->[0]};
        foreach $line (keys %{$self->{'seq'}}) {
            $line = $self->{'seq'}->{$line};
            if (length $line != $off) {
                $self->die("unequal line lengths (expect $off, got @{[length $line]})\n");
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
    return $s;
}


###########################################################################
1;
