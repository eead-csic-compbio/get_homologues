# Copyright (C) 1998-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
# FASTA2 suite ALIGN pairwise method.
###########################################################################
package Bio::Parse::Format::FASTA2::align;

use Bio::Parse::Format::FASTA;
use strict;

use vars qw(
            @ISA

            $ALIGN_START
            $ALIGN_END
);

@ISA = qw(Bio::Parse::Format::FASTA);

#delimit full ALIGN or LALIGN entry
$ALIGN_START             = '^\s*ALIGN calculates';
$ALIGN_END               = $ALIGN_START;

#ALIGN record types
my $ALIGN_ALIGNMENT      = '^(?:\s+\d+\s+|\S+)';    #the ruler or any text
my $ALIGN_ALIGNMENTend   = $ALIGN_END;
my $ALIGN_SUMMARY        = '^\s*\S+\s+identity';
my $ALIGN_SUMMARYend     = $ALIGN_ALIGNMENT;
my $ALIGN_MATCH          = $ALIGN_SUMMARY;
my $ALIGN_MATCHend       = $ALIGN_ALIGNMENTend;
my $ALIGN_HEADER         = $ALIGN_START;
my $ALIGN_HEADERend      = $ALIGN_MATCH;
my $ALIGN_Null           = '^\s*$';#'

#Parse one entry
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #HEADER lines
        if ($line =~ /$ALIGN_HEADER/o) {
            $scan->scan_until($ALIGN_HEADERend);
            $self->push_record('HEADER',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #consume data

        #MATCH lines
        if ($line =~ /$ALIGN_MATCH/o) {
            $scan->scan_until($ALIGN_MATCHend);
            $self->push_record('MATCH',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$ALIGN_Null/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::FASTA2::align::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'program'} = '?';
    $self->{'version'} = '?';
    $self->{'moltype'} = '?';

    $self->{'id1'}     = '?';
    $self->{'desc1'}   = '';
    $self->{'length1'} = 0;

    $self->{'id2'}     = '?';
    $self->{'desc2'}   = '';
    $self->{'length2'} = 0;

    #consume Name lines
    while (defined ($line = $scan->read_line)) {
        #warn $line;

        #program information
        if ($line =~ /^\s*ALIGN/x) {
            $self->{'program'} = 'ALIGN';
            next;
        }

        #ALIGN version information
        if ($line =~ /^
            \s*version\s+(\S+)Please
            /xo) {

            $self->test_args(\$line, $1);
            (
             $self->{'version'},
            ) = ($1);
            next;
        }

        #first sequence: id, description, length
        if ($line =~ /^
            \s*(\S+)\s+              #id
            \s*(.*)\s+               #description (empty?)
            \s*(\d+)\s+(aa|nt)\s+vs  #length
            /xo) {

            $self->test_args(\$line, $1, $3, $4);
            (
             $self->{'id1'},
             $self->{'desc1'},
             $self->{'length1'},
             $self->{'moltype'},
             ) = ($1, (defined $2?$2:''), $3, $4);
            next;
        }

        #second sequence: id, description, length
        if ($line =~ /^
            \s*(\S+)\s+           #id
            \s*(.*)\s+            #description (empty?)
            \s*(\d+)\s+(?:aa|nt)  #length
            /xo) {

            $self->test_args(\$line, $1, $3);
            (
             $self->{'id2'},
             $self->{'desc2'},
             $self->{'length2'},
             ) = ($1, (defined $2?$2:''), $3);
            next;
        }

        #ignore any other text
    }
    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> %s\n",   'program', $self->{'program'};
    $s .= sprintf "$x%20s -> %s\n",   'version', $self->{'version'};
    $s .= sprintf "$x%20s -> %s\n",   'moltype', $self->{'moltype'};
    $s .= sprintf "$x%20s -> %s\n",   'id1',     $self->{'id1'};
    $s .= sprintf "$x%20s -> %s\n",   'desc1',   $self->{'desc1'};
    $s .= sprintf "$x%20s -> %s\n",   'length1', $self->{'length1'};
    $s .= sprintf "$x%20s -> %s\n",   'id2',     $self->{'id2'};
    $s .= sprintf "$x%20s -> %s\n",   'desc2',   $self->{'desc2'};
    $s .= sprintf "$x%20s -> %s\n",   'length2', $self->{'length2'};
    return $s;
}


###########################################################################
package Bio::Parse::Format::FASTA2::align::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        if ($line =~ /$ALIGN_SUMMARY/o) {
            $scan->scan_until($ALIGN_SUMMARYend);
            $self->push_record('SUM',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        if ($line =~ /$ALIGN_ALIGNMENT/o) {
            $scan->scan_until($ALIGN_ALIGNMENTend);
            $self->push_record('ALN',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$ALIGN_Null/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}

###########################################################################
package Bio::Parse::Format::FASTA2::align::MATCH::SUM;

use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'identity'} = '?';
    $self->{'score'}    = '?';

    #consume Name lines
    while (defined ($line = $scan->read_line)) {
        #print $line;

        #ALIGN
        if ($line =~ /^
            \s*($RX_Ureal)%\s+identity
            .*
            score:\s+($RX_Sint)
            /xo) {

            $self->{'identity'} = $1;
            $self->{'score'}    = $2;
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$ALIGN_Null/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> %s\n", 'identity', $self->{'identity'};
    $s .= sprintf "$x%20s -> %s\n", 'score',    $self->{'score'};
    return $s;
}

###########################################################################
package Bio::Parse::Format::FASTA2::align::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA2::MATCH::ALN);


###########################################################################
1;
