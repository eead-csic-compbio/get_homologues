# -*- perl -*-
# Copyright (C) 1996-2011 Nigel P. Brown
# $Id: BLAST1.pm,v 1.12 2013/09/09 21:31:04 npb Exp $

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
package NPB::Parse::Format::BLAST1;

use NPB::Parse::Format::BLAST;
use NPB::Parse::Regexps;

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

@ISA   = qw(NPB::Parse::Format::BLAST);

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

$WARNING_START	  = '^(?:WARNING|NOTE):';
$WARNING_END	  = $NULL;

$WARNINGS_START	  = '^WARNINGS\s+ISSUED:';
$WARNINGS_END	  = $NULL;

$PARAMETERS_START = '^Parameters';
$PARAMETERS_END	  = "^(?:$WARNINGS_START|$ENTRY_START)";

$MATCH_START	  = '^>';
$MATCH_END	  = "^(?:$MATCH_START|$WARNING_START|$PARAMETERS_START|$PARAMETERS_END)";

$RANK_START	  = '^\s+Smallest';
$RANK_MATCH	  = "^(?:$NULL|\\s+Smallest|\\s+Sum|\\s+High|\\s+Reading|\\s*Sequences|[^>].*$RX_Uint\\s+$RX_Ureal\\s+$RX_Uint|$NPB::Parse::Format::BLAST::GCG_JUNK)";
$RANK_END	  = "^(?:$WARNING_START|$PARAMETERS_START|$PARAMETERS_END)";

$HISTOGRAM_START  = '^\s+Observed Numbers';
$HISTOGRAM_END	  = "^(?:$RANK_START|$RANK_END)";

$HEADER_START     = $ENTRY_START;
$HEADER_END       = "^(?:$WARNING_START|$HISTOGRAM_START|$HISTOGRAM_END)";

$SCORE_START      = '^ Score';
$SCORE_END        = "^(?:$SCORE_START|$MATCH_START|$PARAMETERS_START)";


sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);

    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    while (defined ($line = $text->next_line)) {

	#Header lines
	if ($line =~ /$HEADER_START/o) {
	    $text->scan_until($HEADER_END, 'HEADER');
	    next;
	}

	#Histogram lines
	if ($line =~ /$HISTOGRAM_START/o) {
	    $text->scan_until($HISTOGRAM_END, 'HISTOGRAM');
	    next;
	}

	#Rank lines: override $RANK_END definition
	if ($line =~ /$RANK_START/o) {       	      
	    $text->scan_while($RANK_MATCH, 'RANK');
	    next;
	}				       	      
	
	#Hit lines
	if ($line =~ /$MATCH_START/o) {
	    $text->scan_until($MATCH_END, 'MATCH');
	    next;			       	      
	}				       	      
	
	#WARNING lines
	if ($line =~ /$WARNING_START/o) {       	      
	    $text->scan_until($WARNING_END, 'WARNING');
	    next;			       	      
	}				       	      

	#Parameter lines
	if ($line =~ /$PARAMETERS_START/o) {       	      
	    $text->scan_until($PARAMETERS_END, 'PARAMETERS');
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
package NPB::Parse::Format::BLAST1::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::HEADER);


###########################################################################
package NPB::Parse::Format::BLAST1::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::RANK);


###########################################################################
package NPB::Parse::Format::BLAST1::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::MATCH);


###########################################################################
package NPB::Parse::Format::BLAST1::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::BLAST1::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::MATCH::ALN);


###########################################################################
1;
