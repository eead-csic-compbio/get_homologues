# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Base classes for NCBI BLAST1, WashU BLAST2 families.
#
# Handles: BLAST 1.4.x, WashU 2.0x
#
# BLAST (pre NCBI version 2) parsing consists of 6 main record types:
#
#   HEADER        the header text
#   WARNING       optional warning messages
#   HISTOGRAM     the optional scores histogram
#   RANK          the list of ordered high scoring hits
#   MATCH {*}     the set of fragments (HSPs) for a given hit
#     SUM           the summary lines for each hit
#     ALN {*}       each aligned fragment: score + alignment
#   PARAMETERS    the trailer
#
###########################################################################
package Bio::Parse::Format::BLAST1;

use Bio::Parse::Format::BLAST;
use Bio::Util::Regexp;

use strict;

use vars qw(@ISA

            @VERSIONS

            $ENTRY_START
            $ENTRY_END

            $WARNING_START
            $WARNING_END

            $WARNINGS_START
            $WARNINGS_END

            $PARAMETERS_START
            $PARAMETERS_END

            $MATCH_START
            $MATCH_END

            $RANK_START
            $RANK_MATCH
            $RANK_END

            $HISTOGRAM_START
            $HISTOGRAM_END

            $HEADER_START
            $HEADER_END

            $SCORE_START
            $SCORE_END
           );

@ISA   = qw(Bio::Parse::Format::BLAST);

@VERSIONS = (
             '1' => [
                     'BLASTP',
                     'BLASTN',
                     'BLASTX',
                     'TBLASTN',
                     'TBLASTX',
                    ],
            );

my $NULL = '^\s*$';

$ENTRY_START = '^(?:'
    . 'BLASTP'
    . '|'
    . 'BLASTN'
    . '|'
    . 'BLASTX'
    . '|'
    . 'TBLASTN'
    . '|'
    . 'TBLASTX'
    . ')';
$ENTRY_END        = '^WARNINGS\s+ISSUED:';

$WARNING_START    = '^(?:WARNING|NOTE):';
$WARNING_END      = $NULL;

$WARNINGS_START   = '^WARNINGS\s+ISSUED:';
$WARNINGS_END     = $NULL;

$PARAMETERS_START = '^Parameters';
$PARAMETERS_END   = "^(?:$WARNINGS_START|$ENTRY_START)";

$MATCH_START      = '^>';
$MATCH_END        = "^(?:$MATCH_START|$WARNING_START|$PARAMETERS_START|$PARAMETERS_END)";

$RANK_START       = '^\s+Smallest';
$RANK_MATCH       = "^(?:$NULL|\\s+Smallest|\\s+Sum|\\s+High|\\s+Reading|\\s*Sequences|[^>].*$RX_Uint\\s+$RX_Ureal\\s+$RX_Uint|$Bio::Parse::Format::BLAST::GCG_JUNK)";
$RANK_END         = "^(?:$WARNING_START|$PARAMETERS_START|$PARAMETERS_END)";

$HISTOGRAM_START  = '^\s+Observed Numbers';
$HISTOGRAM_END    = "^(?:$RANK_START|$RANK_END)";

$HEADER_START     = $ENTRY_START;
$HEADER_END       = "^(?:$WARNING_START|$HISTOGRAM_START|$HISTOGRAM_END)";

$SCORE_START      = '^ Score';
$SCORE_END        = "^(?:$SCORE_START|$MATCH_START|$PARAMETERS_START)";


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

        #Histogram lines
        if ($line =~ /$HISTOGRAM_START/o) {
            $scan->scan_until($HISTOGRAM_END);
            $self->push_record('HISTOGRAM',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Rank lines: override $RANK_END definition
        if ($line =~ /$RANK_START/o) {
            $scan->scan_while($RANK_MATCH);
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

        #WARNING lines
        if ($line =~ /$WARNING_START/o) {
            $scan->scan_until($WARNING_END);
            $self->push_record('WARNING',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Parameter lines
        if ($line =~ /$PARAMETERS_START/o) {
            $scan->scan_until($PARAMETERS_END);
            $self->push_record('PARAMETERS',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #WARNINGS ISSUED line: ignore
        next    if $line =~ /$WARNINGS_START/o;

        #blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

        #default
        $self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::BLAST1::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST1::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::RANK);


###########################################################################
package Bio::Parse::Format::BLAST1::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST1::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST1::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::MATCH::ALN);


###########################################################################
1;
