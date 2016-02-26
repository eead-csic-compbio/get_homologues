# -*- perl -*-
# Copyright (C) 1998-2015 Nigel P. Brown
# $Id: SIM.pm,v 1.3 2015/06/14 17:09:04 npb Exp $

###########################################################################
# Huang Algorithms SIM
###########################################################################
package NPB::Parse::Format::SIM;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full SIM or LSIM entry
my $SIM_START          = '^Match';
my $SIM_END            = $SIM_START;

#SIM record types
my $SIM_ALNSEP         = '^\*{57}';
my $SIM_SIMMALN        = '^(?:\s+\d+\s+|\S+)';    #the ruler or any text
my $SIM_SIMMALNend     = "(?:$SIM_ALNSEP|$SIM_END)";
my $SIM_SUMMARY        = '^\s*Number';
my $SIM_SUMMARYend     = $SIM_SIMMALN;
my $SIM_MATCH          = $SIM_SUMMARY;
my $SIM_MATCHend       = $SIM_ALNSEP;
my $SIM_HEADER         = $SIM_START;
my $SIM_HEADERend      = $SIM_MATCH;
my $SIM_Null           = '^\s*$';#'


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new SIM instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    my $fixlalign = 0;

    while (defined ($line = <$fh>)) {

	#start of entry
 	if ($line =~ /$SIM_START/o and $offset < 0) {
	    #print "STRT:$line";
	    $offset = $fh->tell - length($line);
	    next;
	}

	#consume rest of stream, until next entry, if any
        if ($line =~ /$SIM_END/o) {
	    $fh->seek(-length($line), 1);
	    #print "STOP:$line";
	    last;
        }
	#print "SEEN:$line";
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::SIM(undef, $text, $offset, $bytes);
}
	    
#Parse one entry
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

	#HEADER lines
	if ($line =~ /$SIM_HEADER/o) {
	    $text->scan_until($SIM_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#MATCH lines		       	      
	if ($line =~ /$SIM_MATCH/o) {
	    $text->scan_until($SIM_MATCHend, 'MATCH');
	    next;			       	      
	}				       	      
	
	#skip LALIGN alignment delimiter
	next    if $line =~ /$SIM_ALNSEP/o;

	#blank line or empty record: ignore
	next    if $line =~ /$SIM_Null/o;

	#default
	$self->warn("unknown fiSSeld: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::SIM::HEADER;

use NPB::Parse::Regexps;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

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

    $self->{'program'} = 'SIM';
    $self->{'version'} = '?';
    $self->{'id1'}     = '?';
    $self->{'desc1'}   = '';
    $self->{'length1'} = 0;
    $self->{'id2'}     = '?';
    $self->{'desc2'}   = '';
    $self->{'length2'} = 0;

    $self->{match}     = '?';
    $self->{mismatch}  = '?';
    $self->{gapopen}   = '?';
    $self->{gapextend} = '?';

    #consume Name lines
    while (defined ($line = $text->next_line)) { 
	#print $line;	

	#scoring parameter line
	#Match   Mismatch   Gap-Open Penalty   Gap-Extension Penalty
	next  if $line =~ /^\s*Match/;

	if ($line =~ /^
	    \s+($RX_Sint)               #match
	    \s+($RX_Sint)               #mismatch
	    \s+($RX_Sint)               #gap-open penalty
	    \s+($RX_Sint)               #gap-extension penalty
	    /xo) {
	    
	    $self->test_args($line, $1, $2, $3, $4);
	    (
	     $self->{'match'},
	     $self->{'mismatch'},
	     $self->{'gapopen'},
	     $self->{'gapextend'},
	     ) = ($1, $2, $3, $4);
	    next;
	}

	#first sequence: id, description
	if ($line =~ /^\s*Upper Sequence:\s+(\S+)\s*(.*)/) {
	    ($self->{'id1'},
	     $self->{'desc1'},
	     ) = ($1, (defined $2?$2:''));
	    next;
	}

	#second sequence: id, description
	if ($line =~ /^\s*Lower Sequence:\s+(\S+)\s*(.*)/) {
	    ($self->{'id2'},
	     $self->{'desc2'},
	     ) = ($1, (defined $2?$2:''));
	    next;
	}
	
	#sequence length
	if ($line =~ /^\s*Length:\s+(\S+)/) {
	    if ($self->{'length1'} == 0) {
		$self->{'length1'} = $1;
	    } else {
		$self->{'length2'} = $1;
	    }
	    next;
	}

	#ignore any other text
    }

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",  'program',       $self->{'program'};
    printf "$x%20s -> %s\n",  'version',       $self->{'version'};
    printf "$x%20s -> %s\n",  'id1',           $self->{'id1'};
    printf "$x%20s -> %s\n",  'desc1',         $self->{'desc1'};
    printf "$x%20s -> %s\n",  'length1',       $self->{'length1'};
    printf "$x%20s -> %s\n",  'id2',           $self->{'id2'};
    printf "$x%20s -> %s\n",  'desc2',         $self->{'desc2'};
    printf "$x%20s -> %s\n",  'length2',       $self->{'length2'};
    printf "$x%20s -> %s\n",  'match',         $self->{'match'};
    printf "$x%20s -> %s\n",  'mismatch',      $self->{'mismatch'};
    printf "$x%20s -> %s\n",  'gap opening',   $self->{'gapopen'};
    printf "$x%20s -> %s\n",  'gap extension', $self->{'gapextend'};
}


###########################################################################
package NPB::Parse::Format::SIM::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid argument list (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);

    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    while (defined ($line = $text->next_line)) {

	if ($line =~ /$SIM_SUMMARY/o) {
	    $text->scan_until($SIM_SUMMARYend, 'SUM');
	    next;
	}

	if ($line =~ /$SIM_SIMMALN/o) {
	    $text->scan_until($SIM_SIMMALNend, 'ALN');
	    next;
	}

	#blank line or empty record: ignore
        next    if $line =~ /$SIM_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

###########################################################################
package NPB::Parse::Format::SIM::MATCH::SUM;

use NPB::Parse::Regexps;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

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

    $self->{'identity'} = '?';
    $self->{'score'}    = '?';

    #consume Name lines
    while (defined ($line = $text->next_line)) { 
	#print $line;

	if ($line =~ /^\s*Number/) {
	    next;
	}
	if ($line =~ /^\s*Similarity Score\s*:\s*(\S+)/) {
	    $self->{score} = $1;
	    next;
	}
	if ($line =~ /^\s*Match Percentage\s*:\s+(\S+)%/) {
	    $self->{'identity'} = $1;
	    next;
	}
	if ($line =~ /^\s*Number of Matches\s*:\s*(\S+)/) {
	    next;
	}
	if ($line =~ /^\s*Number of Mismatches\s*:\s*(\S+)/) {
	    next;
	}
	if ($line =~ /^\s*Total Length of Gaps\s*:\s*(\S+)/) {
	    next;
	}
	if ($line =~ /^\s*Begins at \((\d+)\s*,\s*(\d+)\) and Ends at \((\d+)\s*,\s*(\d+)\)/) {
	    next;
	}

	#blank line or empty record: ignore
        next    if $line =~ /$SIM_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",  'identity', $self->{'identity'};
    printf "$x%20s -> %s\n",  'score',    $self->{'score'};
}

###########################################################################
package NPB::Parse::Format::SIM::MATCH::ALN;

use NPB::Parse::Math;
use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

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

    my ($query_start, $query_stop, $sbjct_start, $sbjct_stop) = (0,0,0,0);
    my ($query, $align, $sbjct, $len) = ('','','');
    my ($x, $depth);

    #first record
    my @tmp = ();

    #warn $text->inspect_stream;

    $line = $text->next_line(1);
    
    if ($line =~ /^\s*[0-9]+[.: ]*$/) {
	$tmp[0] = $line;                  #ruler
	$tmp[1] = $text->next_line(1);    #query sequence
	$tmp[2] = $text->next_line(1);    #match pattern
	$tmp[3] = $text->next_line(1);    #sbjct sequence
    } else {
	$self->die("unexpected line: $line\n");
    }
    
    #initialise query and match start
    $tmp[1] =~ /^\s*(\S+)\s*/;
    $query_start = $1;
    $tmp[3] =~ /^\s*(\S+)\s*/;
    $sbjct_start = $1;

    #determine depth of sequence labels at left: take the shorter
    #since either sequence may begin with a gap
    $tmp[1] =~ /^(\s*\S+\s*)/;
    $depth = length $1;
    $tmp[3] =~ /^(\s*\S+\s*)/;
    $x = length $1;
    $depth = ($depth < $x ? $depth : $x);
    
    #strip width of leading number
    $tmp[0] = substr($tmp[0], $depth);
    $tmp[1] = substr($tmp[1], $depth);
    $tmp[2] = substr($tmp[2], $depth);
    $tmp[3] = substr($tmp[3], $depth);
    
    #warn "#0##$tmp[0]\n";
    #warn "#1##$tmp[1]\n";
    #warn "#2##$tmp[2]\n";
    #warn "#3##$tmp[3]\n";
    
    #query/sbjct/match lines
    #warn "|$tmp[1]|\n|$tmp[2]|\n|$tmp[3]|\n";
    $len = max(length($tmp[1]), length($tmp[3]));
    $tmp[1] .= ' ' x ($len-length $tmp[1])    if length $tmp[1] < $len;
    $tmp[2] .= ' ' x ($len-length $tmp[2])    if length $tmp[2] < $len;
    $tmp[3] .= ' ' x ($len-length $tmp[3])    if length $tmp[3] < $len;

    #first query and sbjct stop
    ($x = $tmp[1]) =~ s/[ -]//g;
    $query_stop = $query_start + length($x) - 1;
    ($x = $tmp[3]) =~ s/[ -]//g;
    $sbjct_stop = $sbjct_start + length($x) - 1;

    #assemble sequences
    $query = $tmp[1];
    $align = $tmp[2];
    $sbjct = $tmp[3];
    
    #remaining records
    while (defined ($line = $text->next_line(1))) {

	next if $line =~ /^\s*$/;

 	@tmp = ();

	if ($line =~ /^\s*[0-9]+[.: ]*$/) {
	    $tmp[0] = $line;                  #ruler
	    $tmp[1] = $text->next_line(1);    #query sequence
	    $tmp[2] = $text->next_line(1);    #match pattern
	    $tmp[3] = $text->next_line(1);    #sbjct sequence
	} else {
	    $self->die("unexpected line: $line\n");
	}
    
	#if ($tmp[1] =~ /^\s*(\d+)/) {
	#    warn "query_stop: old=$query_stop new=$1\n";
	#}
	#if ($tmp[3] =~ /^\s*(\d+)/) {
	#    warn "sbjct_stop: old=$query_stop new=$1\n";
	#}

	#strip width of leading number
	$tmp[0] = substr($tmp[0], $depth);
	$tmp[1] = substr($tmp[1], $depth);
	$tmp[2] = substr($tmp[2], $depth);
	$tmp[3] = substr($tmp[3], $depth);
	
	#warn "#0##$tmp[0]\n";
	#warn "#1##$tmp[1]\n";
	#warn "#2##$tmp[2]\n";
	#warn "#3##$tmp[3]\n";

	#next query and sbjct stop
	($x = $tmp[1]) =~ s/[ -]//g;
	$query_stop = $query_stop + length($x);
	($x = $tmp[3]) =~ s/[ -]//g;
	$sbjct_stop = $sbjct_stop + length($x);

	#assemble sequences
	$query .= $tmp[1];
	$align .= $tmp[2];
	$sbjct .= $tmp[3];
    }
    
    #warn $text->inspect_stream;
    #warn "EXIT ($query_start, $query_stop) ($sbjct_start, $sbjct_stop)\n";

    $self->{'query'} 	   = $query;
    $self->{'align'} 	   = $align;
    $self->{'sbjct'} 	   = $sbjct;

    $self->{'query_start'} = $query_start;
    $self->{'query_stop'}  = $query_stop;

    $self->{'sbjct_start'} = $sbjct_start;
    $self->{'sbjct_stop'}  = $sbjct_stop;

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n",  'query',       $self->{'query'};
    printf "$x%20s -> '%s'\n",  'align',       $self->{'align'};
    printf "$x%20s -> '%s'\n",  'sbjct',       $self->{'sbjct'};
    printf "$x%20s -> %s\n",    'query_start', $self->{'query_start'};
    printf "$x%20s -> %s\n",    'query_stop',  $self->{'query_stop'};
    printf "$x%20s -> %s\n",    'sbjct_start', $self->{'sbjct_start'};
    printf "$x%20s -> %s\n",    'sbjct_stop' , $self->{'sbjct_stop'};
}

###########################################################################
1;
