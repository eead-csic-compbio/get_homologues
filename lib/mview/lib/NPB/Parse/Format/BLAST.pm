# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: BLAST.pm,v 1.21 2015/06/14 17:09:04 npb Exp $

###########################################################################
#
# Base classes for NCBI BLAST1, WashU BLAST2, NCBI BLAST2 families.
#
###########################################################################
package NPB::Parse::Format::BLAST;

use vars qw(@ISA $GCG_JUNK);
use strict;

@ISA = qw(NPB::Parse::Record);

BEGIN { $GCG_JUNK = '^(?:\.\.|\\\\)' }

use NPB::Parse::Format::BLAST1;
use NPB::Parse::Format::BLAST2;
use NPB::Parse::Format::BLAST2_OF7;

my %VERSIONS = (
		@NPB::Parse::Format::BLAST1::VERSIONS,
		@NPB::Parse::Format::BLAST2::VERSIONS,
		@NPB::Parse::Format::BLAST2_OF7::VERSIONS,
	       );

my $NULL = '^\s*$';;

my $ENTRY_START = "(?:"
    . $NPB::Parse::Format::BLAST1::ENTRY_START
    . "|"
    . $NPB::Parse::Format::BLAST2::ENTRY_START
    . "|"
    . $NPB::Parse::Format::BLAST2_OF7::ENTRY_START
    . ")";
my $ENTRY_END   = "(?:"
    . $NPB::Parse::Format::BLAST1::ENTRY_END
    . "|"
    . $NPB::Parse::Format::BLAST2::ENTRY_END
    . "|"
    . $NPB::Parse::Format::BLAST2_OF7::ENTRY_END
    . ")";

my $HEADER_START = "(?:"
    . $NPB::Parse::Format::BLAST1::HEADER_START
    . "|"
    . $NPB::Parse::Format::BLAST2::HEADER_START
    . "|"
    . $NPB::Parse::Format::BLAST2_OF7::HEADER_START
    . ")";
my $HEADER_END   = "(?:"
    . $NPB::Parse::Format::BLAST1::HEADER_END
    . "|"
    . $NPB::Parse::Format::BLAST2::HEADER_END
    . "|"
    . $NPB::Parse::Format::BLAST2_OF7::HEADER_END
    . ")";

my $PARAMETERS_START = "(?:"
    . $NPB::Parse::Format::BLAST1::PARAMETERS_START
    . "|"
    . $NPB::Parse::Format::BLAST2::PARAMETERS_START
    . ")";
my $PARAMETERS_END = "(?:"
    . $NPB::Parse::Format::BLAST1::PARAMETERS_END
    . "|"
    . $NPB::Parse::Format::BLAST2::PARAMETERS_END
    . ")";

my $WARNINGS_START = "(?:"
    . $NPB::Parse::Format::BLAST1::WARNINGS_START
    . "|"
    . $NPB::Parse::Format::BLAST2::WARNINGS_START
    . ")";
my $WARNINGS_END = "(?:"
    . $NPB::Parse::Format::BLAST1::WARNINGS_END
    . "|"
    . $NPB::Parse::Format::BLAST2::WARNINGS_END
    . ")";

my $HISTOGRAM_START = "(?:"
    . $NPB::Parse::Format::BLAST1::HISTOGRAM_START
    . ")";
my $HISTOGRAM_END = "(?:"
    . $NPB::Parse::Format::BLAST1::HISTOGRAM_END
    . ")";

my $RANK_START = "(?:"
    . $NPB::Parse::Format::BLAST1::RANK_START
    . "|"
    . $NPB::Parse::Format::BLAST2::RANK_START
    . ")";
my $RANK_END = "(?:"
    . $NPB::Parse::Format::BLAST1::RANK_END
    . "|"
    . $NPB::Parse::Format::BLAST2::RANK_END
    . ")";

my $MATCH_START = "(?:"
    . $NPB::Parse::Format::BLAST1::MATCH_START
    . "|"
    . $NPB::Parse::Format::BLAST2::MATCH_START
    . ")";
my $MATCH_END = "(?:"
    . $NPB::Parse::Format::BLAST1::MATCH_END
    . "|"
    . $NPB::Parse::Format::BLAST2::MATCH_END
    . ")";

my $WARNING_START = "(?:"
    . $NPB::Parse::Format::BLAST1::WARNING_START
    . "|"
    . $NPB::Parse::Format::BLAST2::WARNING_START
    . ")";
my $WARNING_END = "(?:"
    . $NPB::Parse::Format::BLAST1::WARNING_END
    . "|"
    . $NPB::Parse::Format::BLAST2::WARNING_END
    . ")";

my $SEARCH_START = "(?:"
    . $NPB::Parse::Format::BLAST2::SEARCH_START
    . ")";
my $SEARCH_END = "(?:"
    . $NPB::Parse::Format::BLAST2::SEARCH_END
    . ")";

my $SCORE_START = "(?:"
    . $NPB::Parse::Format::BLAST1::SCORE_START
    . "|"
    . $NPB::Parse::Format::BLAST2::SCORE_START
    . ")";
my $SCORE_END = "(?:"
    . $NPB::Parse::Format::BLAST1::SCORE_END
    . "|"
    . $NPB::Parse::Format::BLAST2::SCORE_END
    . ")";

#Generic get_entry() and new() constructors for all BLAST style parsers:
#determine program and version and coerce appropriate subclass.

#Consume one entry-worth of input on text stream associated with $file and
#return a new BLAST instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my ($type, $prog, $version) = ('NPB::Parse::Format::BLAST', undef, undef);
    my ($savefmt, $saveformat, $saveversion) = ('', '', '');

    while (defined ($line = $parent->{'text'}->getline)) {

        #start of entry
        if ($line =~ /$ENTRY_START/o and $offset < 0) {
            $offset = $parent->{'text'}->startofline;
            #fall through for version tests
        }

	#end of entry
	last  if $line =~ /$ENTRY_END/o;

	#escape iteration if we've found the BLAST type previously
	next  if defined $prog and defined $version;
	
	if ($line =~ /$HEADER_START\s+(\d+)/o) {

	    #read major version
	    $version = $1;

	    #reassess version
	    if ($line =~ /WashU/) {
		#WashU series is almost identical to NCBI BLAST1
		$version = 1;
	    }

	    #read program name and determine output format
	    $line =~ /^(\#)?\s*(\S+)/;
            $savefmt = defined $1 ? '_OF7' : '';
	    $prog = $2;

	    $saveformat  = uc $prog;
	    $saveversion = $version;
		
	    next;
	}
	
    }
    return 0  if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    unless (defined $prog and defined $version) {
	die "get_entry() top-level BLAST parser could not determine program/version\n";
    }

    unless (exists $VERSIONS{$version} and 
	    grep(/^$prog$/i, @{$VERSIONS{$version}}) > 0) {
	die "get_entry() parser for program '$prog' version '$version' not implemented\n";
    }

    if ($prog eq 'PSIBLAST') { # BLAST+
        $prog = 'blastp';
        $saveformat = 'blastp';
    }

    $prog = lc $prog;
    $version =~ s/-/_/g;
    $type = "NPB::Parse::Format::BLAST${version}${savefmt}::$prog";

    #warn "prog= $prog  version= $version  type= $type\n";

    ($prog = $type) =~ s/::/\//g; require "$prog.pm";

    #package $type defines this constructor and coerces to $type

    my $self = $type->new(undef, $parent->{'text'}, $offset, $bytes);

    $self->{'format'}  = $saveformat;
    $self->{'version'} = $saveversion;

    $self;
}

sub new { die "$_[0]::new() virtual function called\n" }


###########################################################################
package NPB::Parse::Format::BLAST::HEADER;

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

    $self->{'full_version'} = '';
    $self->{'version'}      = '',
    $self->{'query'}        = '';
    $self->{'summary'}      = '';

    while ($line = $text->next_line) {

	#blast version info
	if ($line =~ /($HEADER_START\s+(\S+).*)/o) {
	
	    $self->test_args($line, $1, $2);
	    
	    (
	     $self->{'full_version'},
	     $self->{'version'},
	    ) = ($1, $2);
	    
	    next;
	} 

	#query line
	if ($line =~ /^Query=\s+(\S+)?\s*[,;]?\s*(.*)/o) {

	    #no test - either field may be missing
	    $self->{'query'}   = $1  if defined $1;
	    $self->{'summary'} = $2  if defined $2;

	    #strip leading unix path stuff
            $self->{'query'} =~ s/.*\/([^\/]+)$/$1/;

	    next;
	} 

	#ignore any other text
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n", 'version',      $self->{'version'};
    printf "$x%20s -> %s\n", 'full_version', $self->{'full_version'};
    printf "$x%20s -> %s\n", 'query',        $self->{'query'};
    printf "$x%20s -> %s\n", 'summary',      $self->{'summary'};
}


###########################################################################
package NPB::Parse::Format::BLAST::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new { die "$_[0]::new() virtual function called\n" }

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'header', $self->{'header'};
    foreach my $hit (@{$self->{'hit'}}) {
	foreach my $field (sort keys %$hit) {
	    printf "$x%20s -> %s\n", $field,  $hit->{$field};
	}
    }
}


###########################################################################
package NPB::Parse::Format::BLAST::MATCH;

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

    #strand orientations, filled by MATCH::ALN object
    $self->{'orient'} = {};

    while (defined ($line = $text->next_line)) {

	#identifier lines
	if ($line =~ /$MATCH_START/o) {
	    $text->scan_until($SCORE_START, 'SUM');
	    next;
	}

	#scored alignments
	if ($line =~ /$SCORE_START/o) {
	    $text->scan_until($SCORE_END, 'ALN');
	    next;
	}
	
	#blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    foreach my $i (sort keys %{$self->{'orient'}}) {
	printf("$x%20s -> %s\n",
	       "orient $i", scalar @{$self->{'orient'}->{$i}});
    }
}


###########################################################################
package NPB::Parse::Format::BLAST::MATCH::SUM;

use vars qw(@ISA);
use NPB::Parse::Regexps;

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

    $line = $text->scan_lines(0);

    if ($line =~ /^\s*
	>?\s*
	([^\s]+)                            #id
	\s+
	(.*)                                #description
	\s*
	Length\s*=\s*($RX_Uint)             #length
	/xso) {
	
	$self->test_args($line, $1, $3);    #ignore $2

	(
	 $self->{'id'},
	 $self->{'desc'},
	 $self->{'length'},
	) = (NPB::Parse::Record::clean_identifier($1),
	     NPB::Parse::Record::strip_english_newlines($2), $3);
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",   'id',     $self->{'id'};
    printf "$x%20s -> '%s'\n", 'desc',   $self->{'desc'};
    printf "$x%20s -> %s\n",   'length', $self->{'length'};
}


###########################################################################
package NPB::Parse::Format::BLAST::MATCH::ALN;

use vars qw(@ISA);
use NPB::Parse::Math;

@ISA = qw(NPB::Parse::Record);

sub new { die "$_[0]::new() virtual function called\n" }

sub parse_alignment {
    my ($self, $text) = @_;
    my ($query_start, $query_stop, $sbjct_start, $sbjct_stop) = (0,0,0,0);
    my ($query, $align, $sbjct) = ('','','');
    my ($line, $depth, $len);

    my @tmp = ();

    #alignment lines
    while (defined ($line = $text->next_line(1))) {
	
	#blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

	#ignore this line (BLASTN)
	next    if $line =~ /^\s*(?:Plus|Minus) Strand/;

	#first compute sequence indent depth
	if ($line =~ /^(\s*Query\:?\s+(\d+)\s*)/) {
	    $depth = length($1);
	}

	@tmp = ();

	#Query line
	if ($line =~ /^\s*
	    Query\:?
	    \s+
	    (\d+)		#start
	    \s*
	    ([^\d\s]+)?         #sequence
	    \s+
	    (\d+)?		#stop
	    /xo) {

	    #$self->test_args($line, $1, $2, $3);
	    $self->test_args($line, $1);

	    if (defined $2 and defined $3) {
		#normal
		if ($query_start) {
		    $query_stop = $3;
		} else {
		    ($query_start, $query_stop) = ($1, $3);
		}
		$tmp[0] = $2;

	    } else {
		if (! defined $3) {
		    #work around bug in psiblast/BLASTP 2.2.9 [May-01-2004]
		    #which introduces white space into the alignment,
		    #sometimes with no sequence at all, and sometimes no
		    #trailing position.
		    ($query_start, $query_stop) = ($1, $1);
		}
		if (! defined $2) {
		    $tmp[0] = '';
		} else {
		    $tmp[0] = $2;
		}
	    }
	    #fall through

	} else {
	    $self->warn("expecting 'Query' line: [$line]");
	    next;
	}
	
	#force read of match line, but note:
	#PHI-BLAST has an extra line - ignore for the time being
	$line = $text->next_line(1);
	$line = $text->next_line(1)  if $line =~ /^pattern/;

	#alignment line
	$tmp[1] = '';
	$tmp[1] = substr($line, $depth)  if length($line) > $depth;

	#force read of Sbjct line
	$line = $text->next_line;

	#Sbjct line
	if ($line =~ /^\s*
	    Sbjct\:?
	    \s+
	    (\d+)               #start
	    \s*
	    ([^\d\s]+)?		#sequence
	    \s+
	    (\d+)?	        #stop
	    /xo) {
	    
	    #$self->test_args($line, $1, $2, $3);
	    $self->test_args($line, $1);

	    if (defined $2 and defined $3) {
		#normal
		if ($sbjct_start) {
		    $sbjct_stop = $3;
		} else {
		    ($sbjct_start, $sbjct_stop) = ($1, $3);
		}
		$tmp[2] = $2;

	    } else {
		if (! defined $3) {
		    #work around bug in psiblast/BLASTP 2.2.9 [May-01-2004]
		    #which introduces white space into the alignment,
		    #sometimes with no sequence at all, and sometimes no
		    #trailing position.
		    ($sbjct_start, $sbjct_stop) = ($1, $1);
		}
		if (! defined $2) {
		    $tmp[2] = '';
		} else {
		    $tmp[2] = $2;
		}
	    }
	    #fall through

	} else {
	    $self->warn("expecting 'Sbjct' line: [$line]");
	    next;
	}

	#query/match/sbjct lines
	#warn "|$tmp[0]|\n|$tmp[1]|\n|$tmp[2]|\n";
	$len = max(length($tmp[0]), length($tmp[2]));
	if (length $tmp[0] < $len) {
	    $tmp[0] .= ' ' x ($len-length $tmp[0]);
	}
	if (length $tmp[1] < $len) {
	    $tmp[1] .= ' ' x ($len-length $tmp[1]);
	}
	if (length $tmp[2] < $len) {
	    $tmp[2] .= ' ' x ($len-length $tmp[2]);
	}

	#work around bug in psiblast/BLASTP 2.2.9 [May-01-2004] with
	#white space sequence padding problem.
	my $qlen = length($tmp[0]);
	my $alen = length($tmp[1]);
	if ($alen > $qlen) {
	    #warn "before[$tmp[1]]\n";
	    $tmp[1] = substr($tmp[1], 0, $qlen);
	    #warn "after [$tmp[1]]\n";
	}
	
	#append sequence strings
	$query .= $tmp[0];
	$align .= $tmp[1];
	$sbjct .= $tmp[2];
    }

    if (length($query) != length($align) or length($query) != length($sbjct)) {
	#warn "Q: ", length($query), "\n";
	#warn "A: ", length($align), "\n";
	#warn "S: ", length($sbjct), "\n";
	$self->warn("unequal Query/Align/Sbjct lengths:\n$query\n$align\n$sbjct\n");
    }

    $self->{'query'} 	    = $query;
    $self->{'align'} 	    = $align;
    $self->{'sbjct'} 	    = $sbjct;
    $self->{'query_start'}  = $query_start;
    $self->{'query_stop'}   = $query_stop;
    $self->{'sbjct_start'}  = $sbjct_start;
    $self->{'sbjct_stop'}   = $sbjct_stop;

    #use sequence numbering to get orientations
    $self->{'query_orient'} =
        $self->{'query_start'} > $self->{'query_stop'} ? '-' : '+';
    $self->{'sbjct_orient'} =
        $self->{'sbjct_start'} > $self->{'sbjct_stop'} ? '-' : '+';

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'query',        $self->{'query'};
    printf "$x%20s -> '%s'\n", 'align',        $self->{'align'};
    printf "$x%20s -> '%s'\n", 'sbjct',        $self->{'sbjct'};
    printf "$x%20s -> %s\n",   'query_orient', $self->{'query_orient'};
    printf "$x%20s -> %s\n",   'query_start',  $self->{'query_start'};
    printf "$x%20s -> %s\n",   'query_stop',   $self->{'query_stop'};
    printf "$x%20s -> %s\n",   'sbjct_orient', $self->{'sbjct_orient'};
    printf "$x%20s -> %s\n",   'sbjct_start',  $self->{'sbjct_start'};
    printf "$x%20s -> %s\n",   'sbjct_stop',   $self->{'sbjct_stop'};
}


###########################################################################
package NPB::Parse::Format::BLAST::WARNING;

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
    $text = new NPB::Parse::Record_Stream($self, 0, undef);

    $line = $text->scan_lines(0);

    $self->{'warning'} = NPB::Parse::Record::strip_english_newlines($line);

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'warning', $self->{'warning'};
}


###########################################################################
package NPB::Parse::Format::BLAST::HISTOGRAM;

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
    $text = new NPB::Parse::Record_Stream($self, 0, undef);

    $line = $text->scan_lines(0);

    ($self->{'histogram'} = $line) =~ s/\s*$/\n/;

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'histogram', $self->{'histogram'};
}


###########################################################################
package NPB::Parse::Format::BLAST::PARAMETERS;

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
    $text = new NPB::Parse::Record_Stream($self, 0, undef);

    $line = $text->scan_lines(0);

    ($self->{'parameters'} = $line) =~ s/\s*$/\n/;

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'parameters', $self->{'parameters'};
}


###########################################################################
1;
