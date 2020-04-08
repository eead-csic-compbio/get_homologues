# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Handles: fasta,fastx,fasty,tfasta,tfastx,tfasty,tfastxy   3.x
#
###########################################################################
package Bio::Parse::Format::FASTA3;

use Bio::Parse::Format::FASTA;
use strict;

use vars qw(
            @ISA

            @VERSIONS

            $NULL

            $ENTRY_START
            $ENTRY_END

            $HEADER_START
            $HEADER_END

            $RANK_START
            $RANK_END

            $TRAILER_START
            $TRAILER_END

            $MATCH_START
            $MATCH_END

            $SUM_START
            $SUM_END

            $ALN_START
            $ALN_END
);

@ISA   = qw(Bio::Parse::Format::FASTA);

@VERSIONS = (
             '3' => [
                     'FASTA',
                     'FASTX',
                     'FASTY',
                     'TFASTA',
                     'TFASTX',
                     'TFASTY',
                     'TFASTXY',
                    ],
            );

$NULL  = '^\s*$';

$ENTRY_START   = '(?:^\s*'
    . 'FASTA +searches +a +protein +or +DNA +sequence +data +bank'
    . '|'
    . 'FAST[XY] +compares +a +DNA +sequence +to +a +protein +sequence +data +bank'
    . '|'
    . 'TFAST[AXY] +compares +a +protein +to +a +translated +DNA +data +bank'
    . '|'
    . 'TFASTXY +compares +a +protein +to +a +translated +DNA +data +bank'
    . '|'
    . '\S+\s*[,:]\s+\d+\s+(?:aa|nt)'
    . ')';
$ENTRY_END     = 'Function used was';

$HEADER_START  = $ENTRY_START;
$HEADER_END    = '^The\s+best\s+scores\s+are:';  #3.2t05 added a space!

$RANK_START    = $HEADER_END;
$RANK_END      = $NULL;

$TRAILER_START = '^\d+\s+residues\s+in\s+\d+\s+query';
$TRAILER_END   = $ENTRY_END;

$MATCH_START   = '^>>\S';
$MATCH_END     = "(?:$MATCH_START|$TRAILER_START|$ENTRY_END)";

$SUM_START     = $MATCH_START;
$SUM_END       = $NULL;

$ALN_START     = '^(?:\s+\d+\s+|\s+$)';  #the ruler
$ALN_END       = $MATCH_END;

#Parse one entry: generic for all FASTA3
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #Header lines
        if ($line =~ /$HEADER_START/o) {
            $scan->scan_until($HEADER_END);
            $self->push_record('HEADER',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Rank lines
        if ($line =~ /$RANK_START/o) {
            $scan->scan_until($RANK_END);
            $self->push_record('RANK',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Hit lines
        if ($line =~ /$MATCH_START/o) {
            $scan->scan_until($MATCH_END);
            $self->push_record('MATCH',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Trailer lines
        if ($line =~ /$TRAILER_START/o) {
            $scan->scan_until_inclusive($TRAILER_END);
            $self->push_record('TRAILER',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #end of FASTA job
        next    if $line =~ /$ENTRY_END/o;

        #blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::FASTA3::HEADER;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::HEADER);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'query'}     = '';
    $self->{'queryfile'} = '';

    while (defined ($line = $scan->read_line)) {

        if ($line =~ /^\s*(version\s+(\S+).*)/) {
            $self->{'full_version'} = $1;
            $self->{'version'}      = $2;
            next;
        }

        #fasta < 3.4
        if ($line =~ /^\s*([^\s>]+)\s+(\d+)\s+(?:aa|nt)/o) {
            $self->test_args(\$line, $1, $2);
            (
             $self->{'queryfile'},
             $self->{'length'},
            ) = ($1, $2);
            $self->{'queryfile'} =~ s/,$//;
            next;
        }

        #fasta 3.4
        if ($line =~ /^\s*Query library\s+(\S+)/o) {
            $self->test_args(\$line, $1);
            $self->{'queryfile'} = $1;
            $self->{'queryfile'} =~ s/,$//;
            next;
        }
        if ($line =~ /^\s*\d+>+(\S+).*-\s+(\d+)\s+(?:aa|nt)/) {
            $self->test_args(\$line, $1, $2);
            (
             $self->{'query'},
             $self->{'length'},
            ) = (clean_identifier($1), $2);
            next;
        }

        if ($line =~ /^(\d+)\s+residues\s+in\s+(\d+)\s+sequences/) {
            $self->test_args(\$line, $1,$2);
            (
             $self->{'residues'},
             $self->{'sequences'},
            ) = ($1, $2);
            next;
        }

        #ignore any other text
    }

    if (! defined $self->{'full_version'} ) {
        #can't determine version: hardwire one!
        $self->{'full_version'} = 'looks like FASTA 3';
        $self->{'version'}      = '3.4';
        $self->warn("unknown FASTA series 3 version - using $self->{'version'}");
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA3::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH::ALN);


###########################################################################
1;
