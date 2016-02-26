# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: FASTA2.pm,v 1.19 2015/06/14 17:09:04 npb Exp $

###########################################################################
#
# Handles: fasta, tfastx  2.x
#
###########################################################################
package NPB::Parse::Format::FASTA2;

use NPB::Parse::Format::FASTA;
use NPB::Parse::Format::FASTA2::align;
use NPB::Parse::Format::FASTA2::lalign;

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

@ISA   = qw(NPB::Parse::Format::FASTA);

@VERSIONS = ( 
	     '2' => [
		     'FASTA',
		     'TFASTX',
		     'ALIGN',
		     'LALIGN',
		     'SSEARCH',
		    ],
	    );

$NULL  = '^\s*$';#for emacs';

$ENTRY_START   = '(?:'
    . '^\s*FASTA searches a protein or DNA sequence data bank'
    . '|'
    . '^\s*TFASTX translates and searches a DNA sequence data bank'
    . '|'
    . '^\s*\S+\s*[,:]\s+\d+\s+(?:aa|nt)'
    . '|'
    . 'SSEARCH searches a sequence database'
    . '|'
    . $NPB::Parse::Format::FASTA2::align::ALIGN_START
    . '|'
    . $NPB::Parse::Format::FASTA2::lalign::ALIGN_START
    . ')';
$ENTRY_END     = '(?:'
    . 'Library scan:'
    . '|'
    . $NPB::Parse::Format::FASTA2::align::ALIGN_END
    . '|'
    . $NPB::Parse::Format::FASTA2::lalign::ALIGN_END
    . ')';
    
$HEADER_START  = $ENTRY_START;
$HEADER_END    = '^The best scores are:'; 
               
$RANK_START    = $HEADER_END;
$RANK_END      = $NULL;
               
$TRAILER_START = $ENTRY_END;
$TRAILER_END   = $ENTRY_END;

$MATCH_START   = '^>*\S{7}.*\(\d+ (?:aa|nt)\)';
$MATCH_END     = "(?:$MATCH_START|$ENTRY_END)";

$SUM_START     = $MATCH_START;
$SUM_END       = $NULL;
       
$ALN_START     = '^(?:\s+\d+\s+|\s+$)';  #the ruler
$ALN_END       = $MATCH_END;

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA2::HEADER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::HEADER);

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

    $self->{'query'}     = '';
    $self->{'queryfile'} = '';

    while (defined ($line = $text->next_line)) {
    
	if ($line =~ /^\s*(version\s+(\S+).*)/) {
	    $self->{'full_version'} = $1;
	    $self->{'version'}      = $2;
	    next;
	}

	if ($line =~ /^\s*v((\d+(?:\.\d+)?)\s*\S+.*)/) {
	    $self->{'full_version'} = $1;
	    $self->{'version'}      = $2;
	    next;
	}

	if ($line =~ /^\s*>?(\S+).*\:\s+(\d+)\s+(?:aa|nt)/) {
	    if ($self->{'queryfile'} eq '') {
		$self->{'queryfile'} = $1;
		$self->{'queryfile'} =~ s/,$//;
	    } else {
		$self->test_args($line, $1);
		$self->{'query'} = NPB::Parse::Record::clean_identifier($1);
	    }
	    $self->{'length'} = $2;
	    next;
	}

	if ($line =~ /^(\d+)\s+residues\s+in\s+(\d+)\s+sequences/) {
	    $self->test_args($line, $1,$2);
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
	$self->{'full_version'} = 'looks like FASTA 2';
	$self->{'version'}      = '2.1';
        $self->warn("unknown FASTA version - using $self->{'version'}");    }

    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA2::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA1::MATCH::ALN);


###########################################################################
1;
