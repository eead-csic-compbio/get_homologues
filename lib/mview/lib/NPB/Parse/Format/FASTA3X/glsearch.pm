# -*- perl -*-
# Copyright (C) 2013-2015 Nigel P. Brown
# $Id: glsearch.pm,v 1.2 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::FASTA3X::glsearch;

use NPB::Parse::Format::FASTA3X;

use strict;
use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3X);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA3X::glsearch::HEADER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA3X::HEADER);

sub new {
    my $self = shift;
    $self = $self->SUPER::new(@_);
    #assume the query identifier is the same as the query filename
    if ($self->{query} eq '' and $self->{queryfile} ne '') {
	$self->{query} = $self->{queryfile};
    }
    return $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3X::glsearch::RANK;

use NPB::Parse::Regexps;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::RANK);

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

    #ranked search hits
    while (defined ($line = $text->next_line)) {
	
	next    if $line =~ /$NPB::Parse::Format::FASTA3X::RANK_START/o;

	if($line =~ /^
	   \s*
	   ([^\s]+)                #id
	   \s+
	   (.*)                    #description: possibly empty
	   \s+
	   \(\s*(\S+)\)            #hit sequence length
	   \s*
	   (?:\[(\S)\])?           #frame
	   \s*
	   ($RX_Sint)              #n-w score
	   \s+
	   (\S+)                   #bits
	   \s+
	   (\S+)                   #E-value
	   \s*
	   $/xo) {
	    
	    $self->test_args($line, $1, $3, $5,$6,$7);
	    
	    push(@{$self->{'hit'}},
		 { 
		  'id'      => NPB::Parse::Record::clean_identifier($1),
		  'desc'    => $2,
		  #ignore $3
		  'frame'   => NPB::Parse::Format::FASTA::parse_frame($4),
		  'orient'  => NPB::Parse::Format::FASTA::parse_orient($4),
		  'score'   => $5,
		  'bits'    => $6,
		  'expect'  => $7,
		 });
	    next;
	}
    
	#blank line or empty record: ignore
	next    if $line =~ /$NPB::Parse::Format::FASTA3X::NULL/o;
	
	#default
	$self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3X::glsearch::TRAILER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::TRAILER);


###########################################################################
package NPB::Parse::Format::FASTA3X::glsearch::MATCH;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::MATCH);


###########################################################################
package NPB::Parse::Format::FASTA3X::glsearch::MATCH::SUM;

use NPB::Parse::Regexps;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA::MATCH::SUM);

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

    my $lines = $text->scan_until_inclusive('^\s*global/local');

    if ($lines =~ /^
	>*
	(\S+)                      	  #id
	\s+
	(.*)                       	  #description: possibly empty
	\s+
	\(\s*(\d+)\s*(?:aa|nt)\)   	  #length
	\s*
        (rev-comp)?                       #frame
        \s+
        n-w\s+opt:\s*($RX_Sint)       	  #n-w opt
        \s+
        Z-score:\s*(\S+)           	  #Z-score
        \s+
        bits:\s*(\S+)              	  #bits
        \s+
        E\((?:\d+)?\):\s*(\S+)         	  #expect
        \s+
        global\/local\s+score:\s*($RX_Sint);   #n-w score (same as opt?)
        \s+
        (\S+)%\s+identity                 #percent identity
        \s+
        \((\S+)%\s+(?:similar|ungapped)\) #percent similarity
        \s+in\s+
        (\d+)\s*(?:aa|nt)\s+overlap       #alignment overlap
        \s+
        \((\S+)\)                         #alignment ranges
        \s*
	$/xso) {

	$self->test_args($line, $1, $3, $5,$6,$7,$8,$9,$10,$11,$12,$13);

	(
	 $self->{'id'},
	 $self->{'desc'},
	 $self->{'length'},
	 $self->{'frame'},
	 $self->{'orient'},
	 $self->{'opt'},
	 $self->{'zscore'},
	 $self->{'bits'},
	 $self->{'expect'},
	 $self->{'score'},
	 $self->{'id_percent'},
	 $self->{'sim_percent'},
	 $self->{'overlap'},
	 $self->{'ranges'},
	) = (
	    NPB::Parse::Record::clean_identifier($1),
	    NPB::Parse::Record::strip_english_newlines($2),
	    $3,
	    NPB::Parse::Format::FASTA::parse_frame($4),
	    NPB::Parse::Format::FASTA::parse_orient($4),
	    $5, $6, $7, $8, $9, $10, $11, $12, $13,
	);

    } elsif ($lines =~ /^>--/) {  #alternative alignment
        my $sib = $self->get_sibling(0);
        (
         $self->{'id'},
         $self->{'desc'},
         $self->{'length'},
         $self->{'alternative'},
        ) =
        (
         $sib->{'id'},
         $sib->{'desc'},
         $sib->{'length'},
         1,  #true
        );

    } else {
	$self->warn("unknown field: $lines");
    }

    $self;
}


###########################################################################
package NPB::Parse::Format::FASTA3X::glsearch::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA3X::MATCH::ALN);

sub query_orient { $_[0]->get_summary->{'orient'} }
sub sbjct_orient { '+' }


###########################################################################
1;
