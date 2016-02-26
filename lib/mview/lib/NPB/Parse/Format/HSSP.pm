# -*- perl -*-
# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: HSSP.pm,v 1.26 2015/06/14 17:09:04 npb Exp $

###########################################################################
#
# HSSP parsing consists of 5 levels of record:
#
# HEADER
#    Crudely presents each simple subrecord as a field of the same name,
#    except the HSSP line which is split for version numbers.
#
# PROTEIN
#    The ranking of proteins against the test sequence:
#
#    'ranking' is an array of records, one per line of input ordered by NR
#    with corresponding field names.
#
#      get_record() method returns this array, or just the supplied record
#      numbers counting from 1.
#
# ALIGNMENT
#    The structural information about the test sequence and the alignments:
#
#    'structure' is an array ordered by SeqNo holding the structural data
#    for the test sequences. This data repeats every '## ALIGNMENT' block,
#    so is only read once:
#
#      get_chains() returns a list of chain names.
#
#      get_record() returns a list of records for the structure fields (the
#      alignment portion is not returned), or just the supplied record
#      numbers counting from 1.
#
#      get_query() returns the test sequence string or just that of the
#      specified chain.
#
#      get_dssp() returns the DSSP assignment or just that of the specified
#      chain.
#
#    'alignment' is an array ordered by SeqNo holding the total alignment
#    data over all '## ALIGNMENT' blocks:
#
#      get_sequence() returns the sequence string for the n'th sequence
#      counting from 1 or just that of the specified chain.
#
# PROFILE
#    Subrecord parser not implemented.
#
# INSERTION
#    Subrecord parser not implemented.
#
###########################################################################
package NPB::Parse::Format::HSSP;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full HSSP entry
my $HSSP_START          = '^HSSP';
my $HSSP_END            = '^\/\/';

#HSSP record types
my $HSSP_HEADER         = $HSSP_START;
my $HSSP_HEADERend      = '^## PROTEINS';
my $HSSP_PROTEIN        = '^## PROTEINS';
my $HSSP_PROTEINend     = '^## ALIGNMENTS';
my $HSSP_ALIGNMENT      = '^## ALIGNMENTS';
my $HSSP_ALIGNMENTend   = '^## SEQUENCE';
my $HSSP_PROFILE        = '^## SEQUENCE';
my $HSSP_PROFILEend     = '^## INSERTION';
my $HSSP_INSERTION      = '^## INSERTION';
my $HSSP_INSERTIONend   = $HSSP_END;
my $HSSP_Null           = '^\s*$';#'


#Consume one entry-worth of input on text stream associated with $file and
#return a new HSSP instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {
	
	#start of entry
	if ($line =~ /$HSSP_START/o and $offset < 0) {
	    $offset = $parent->{'text'}->startofline;
	    next;
	}

	#end of entry
	if ($line =~ /$HSSP_END/o) {
	    last;
	}
    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::HSSP(undef, $parent->{'text'}, $offset, $bytes);
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
	if ($line =~ /$HSSP_HEADER/o) {
	    $text->scan_until($HSSP_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#PROTEIN lines	       	      
	if ($line =~ /$HSSP_PROTEIN/o) {       	      
	    $text->scan_until($HSSP_PROTEINend, 'PROTEIN');
	    next;			       	      
	}				       	      
	
	#ALIGNMENT lines		       	      
	if ($line =~ /$HSSP_ALIGNMENT/o) {       	      
	    $text->scan_until($HSSP_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#PROFILE lines		       	      
	if ($line =~ /$HSSP_PROFILE/o) {       	      
	    $text->scan_until($HSSP_PROFILEend, 'PROFILE');
	    next;			       	      
	}				       	      
	
	#INSERTION lines		       	      
	if ($line =~ /$HSSP_INSERTION/o) {       	      
	    $text->scan_until($HSSP_INSERTIONend, 'INSERTION');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	if ($line =~ /$HSSP_Null/o) {
	    next;
	}

	#terminal line: ignore
	if ($line =~ /$HSSP_END/o) {
	    next;
	}

	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::HSSP::HEADER;

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

    my $tmp;

    while (defined ($line = $text->next_line)) {
    
	if ($line =~ /^HSSP\s+(.*VERSION\s+(.*)\s*)/o) {
	    $self->test_args($line, $1, $2);
	    (
	     $self->{'full_version'},
	     $self->{'version'},
	    ) = ($1, $2);
	    chomp $self->{'full_version'};
	    next;
	} 

	if ($line =~ /^PDBID\s+(.*)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'pdbid'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^DATE\s+file generated on\s+(\S+)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'date'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^SEQBASE\s+(.*)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'seqbase'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^PARAMETER/) {
	    $line = $text->scan_while('^PARAMETER');
	    $tmp  = new NPB::Parse::Record_Stream($self, 11, \$line);
	    $line = $tmp->scan_lines(0);
	    chomp $line;
	    $line = NPB::Parse::Record::strip_english_newlines($line);
	    #fix typo in maxhom output
	    $line =~ s/ :/: /g;
	    #add fullstops and trim whitespace
	    $line =~ s/([^:]+):\s*(\S+)\s*/$1: $2. /g;
	    $self->{'parameter'} = $line;
	    next;
	} 

	if ($line =~ /^THRESHOLD\s+(.*)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'threshold'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^REFERENCE\s+(.*)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'reference'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^CONTACT\s+(.*)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'contact'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^AVAILABLE/) {
	    $line = $text->scan_while('^AVAILABLE');
	    $tmp  = new NPB::Parse::Record_Stream($self, 11, \$line);
	    $line = $tmp->scan_lines(0);
	    $self->{'available'} = NPB::Parse::Record::strip_english_newlines($line);
	    next;
	} 

	if ($line =~ /^HEADER/) {
	    $line = $text->scan_while('^HEADER');
	    $tmp  = new NPB::Parse::Record_Stream($self, 11, \$line);
	    $line = $tmp->scan_lines(0);
	    $self->{'header'} = NPB::Parse::Record::strip_english_newlines($line);
	    next;
	} 

	if ($line =~ /^COMPND\s+(.*)/) {
	    $line = $text->scan_while('^COMPND');
	    $tmp  = new NPB::Parse::Record_Stream($self, 11, \$line);
	    $line = $tmp->scan_lines(0);
	    $self->{'compnd'} = NPB::Parse::Record::strip_english_newlines($line);
	    next;
	} 

	if ($line =~ /^SOURCE\s+(.*)/) {
	    $line = $text->scan_while('^SOURCE');
	    $tmp  = new NPB::Parse::Record_Stream($self, 11, \$line);
	    $line = $tmp->scan_lines(0);
	    $self->{'source'} = NPB::Parse::Record::strip_english_newlines($line);
	    next;
	} 

	if ($line =~ /^AUTHOR\s+(.*)/) {
	    $line = $text->scan_while('^AUTHOR');
	    $tmp  = new NPB::Parse::Record_Stream($self, 11, \$line);
	    $line = $tmp->scan_lines(0);
	    $self->{'author'} = NPB::Parse::Record::strip_english_newlines($line);
	    next;
	} 

	if ($line =~ /^SEQLENGTH\s+(\d+)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'seqlength'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^NCHAIN\s+(\d+)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'nchain'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^KCHAIN\s+(\d+)\s+chain\(s\) used here\s*;\s*chain\(s\)\s*:\s*(.*)/) {
	    $self->test_args($line, $1, $2);
	    (
	     $self->{'kchain'},
	     $self->{'chainname'},
	    ) = ($1, [ split(/,\s*/, $2) ]);
	    next;
	} 

	if ($line =~ /^NALIGN\s+(\d+)/) {
	    $self->test_args($line, $1);
	    (
	     $self->{'nalign'},
	    ) = ($1);
	    next;
	} 

	if ($line =~ /^NOTATION\s+(.*)/) {
	    $line = $text->scan_while('^NOTATION');
	    $tmp  = new NPB::Parse::Record_Stream($self, 11, \$line);
	    $line = $tmp->scan_lines(0);
	    chomp $line;
	    $self->{'notation'} = $line;

	    next;
	}

	if ($line =~ /$HSSP_Null/) {
	    next;
	}

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'version',        $self->{'version'};
    printf "$x%20s -> '%s'\n", 'full_version',   $self->{'full_version'};
    printf "$x%20s -> %s\n",   'pdbid',          $self->{'pdbid'};
    printf "$x%20s -> %s\n",   'date',           $self->{'date'};
    printf "$x%20s -> %s\n",   'seqbase',        $self->{'seqbase'};
    printf "$x%20s -> %s\n",   'parameter',      $self->{'parameter'};
    printf "$x%20s -> %s\n",   'threshold',      $self->{'threshold'};
    printf "$x%20s -> %s\n",   'reference',      $self->{'reference'};
    printf "$x%20s -> %s\n",   'contact',        $self->{'contact'};
    printf "$x%20s -> %s\n",   'available',      $self->{'available'};
    printf "$x%20s -> %s\n",   'header',         $self->{'header'};
    printf "$x%20s -> %s\n",   'compnd',         $self->{'compnd'};
    printf "$x%20s -> %s\n",   'source',         $self->{'source'};
    printf "$x%20s -> %s\n",   'author',         $self->{'author'};
    printf "$x%20s -> %s\n",   'seqlength',      $self->{'seqlength'};
    printf "$x%20s -> %s\n",   'nchain',         $self->{'nchain'};
    printf "$x%20s -> %s\n",   'kchain',         $self->{'kchain'};
    printf "$x%20s -> %s\n",   'chainname',     
        "[" . join(',',@{$self->{'chainname'}}) ."]";
    printf "$x%20s -> %s\n",   'nalign',         $self->{'nalign'};
    printf "$x%20s -> %s\n",   'notation',       $self->{'notation'};
}


###########################################################################
package NPB::Parse::Format::HSSP::PROTEIN;

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

    my ($cut, $part1, $part2, $data);

    #(## PROTEINS : EMBL/SWISSPROT ...) line
    $line = $text->next_line;

    #(  NR.    ID         STRID   %IDE ...) line
    $line = $text->next_line;
    $cut   = index($line, '%IDE');

    $self->{'ranking'} = [];

    while (defined ($line = $text->next_line)) {
	
	$part1 = substr($line, 0, $cut+1);
	$part2 = substr($line, $cut);
	
	$data = {};

	if ($part1 =~ /^
	    \s*
	    (\d+)           #NR.
	    \s*:\s*
	    (\S+)           #ID
	    \s*
	    (\S+)?          #STRID
	    /xo) {

	    $self->test_args($line, $1, $2, "$3");

	    (
	     $data->{'nr'},
	     $data->{'id'},
	     $data->{'strid'},
	    ) = ($1, $2, $3 eq 0 ? '' : $3);

	    if ($part2 =~ /^
		\s*
		(\S+)       #%IDE
		\s+
		(\S+)       #%WSIM
		\s+
		(\S+)       #IFIR
		\s+
		(\S+)       #ILAS
		\s+
		(\S+)       #JFIR
		\s+
		(\S+)       #JLAS
		\s+
		(\S+)       #LALI
		\s+
		(\S+)       #NGAP
		\s+
		(\S+)       #LGAP
		\s+
		(\S+)       #LSEQ2
		\s\s
		(.{11})     #ACCNUM
		(.*)?       #PROTEIN
		/xo) {

		$self->test_args($line, $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12);
		(
		 $data->{'%ide'},
		 $data->{'%wsim'},
		 $data->{'ifir'},
		 $data->{'ilas'},
		 $data->{'jfir'},
		 $data->{'jlas'},
		 $data->{'lali'},
		 $data->{'ngap'},
		 $data->{'lgap'},
		 $data->{'lseq2'},
		 $data->{'protein'},
		) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12);

		$data->{'accnum'} = NPB::Parse::Record::strip_trailing_space($11);

	    }

	    push @{$self->{'ranking'}}, $data;
	    
	    next;
	}

	#default
	$self->warn("unknown field: $line");
    }


    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    for (my $i=0; $i<@{$self->{'ranking'}}; $i++) {
	printf "$x%20s -> %s\n", 'nr',      $self->{'ranking'}->[$i]->{'nr'};
	printf "$x%20s -> %s\n", 'id',      $self->{'ranking'}->[$i]->{'id'};
	printf "$x%20s -> %s\n", '%ide',    $self->{'ranking'}->[$i]->{'%ide'};
	printf "$x%20s -> %s\n", '%wsim',   $self->{'ranking'}->[$i]->{'%wsim'};
	printf "$x%20s -> %s\n", 'ifir',    $self->{'ranking'}->[$i]->{'ifir'};
	printf "$x%20s -> %s\n", 'ilas',    $self->{'ranking'}->[$i]->{'ilas'};
	printf "$x%20s -> %s\n", 'jfir',    $self->{'ranking'}->[$i]->{'jfir'};
	printf "$x%20s -> %s\n", 'jlas',    $self->{'ranking'}->[$i]->{'jlas'};
	printf "$x%20s -> %s\n", 'lali',    $self->{'ranking'}->[$i]->{'lali'};
	printf "$x%20s -> %s\n", 'ngap',    $self->{'ranking'}->[$i]->{'ngap'};
	printf "$x%20s -> %s\n", 'lgap',    $self->{'ranking'}->[$i]->{'lgap'};
	printf "$x%20s -> %s\n", 'lseq2',   $self->{'ranking'}->[$i]->{'lseq2'};
	printf "$x%20s -> %s\n", 'protein', $self->{'ranking'}->[$i]->{'protein'};
    }
}

#returns an ordered array of rank records, or just those requested in the
#order requested. record numbers are counted from 1.
sub get_record {
    my $self = shift;
    my $i;
    my @tmp = ();
    if (@_) {
	foreach $i (@_) {
	    if (defined $self->{'ranking'}->[$i]) {
		push @tmp, $self->{'ranking'}->[$i];
	    }
	}
    } else {
	for ($i=0; $i < @{$self->{'ranking'}}; $i++) {
	    push @tmp, $self->{'ranking'}->[$i];
	}
    }
    @tmp;
}


###########################################################################
package NPB::Parse::Format::HSSP::ALIGNMENT;

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

    my ($lo, $hi, $cut, $part2, $data, $chaincount, $chain1, $chain2, $i);

    #construct an array of query structural properties ordered by SeqNo
    $self->{'structure'}  = [];

    #construct an array of alignments ordered by SeqNo
    $self->{'alignment'}  = [];

    #hash of indices into 'structure' array, keyed by chainname
    $self->{'chain_hash'} = {};

    #ordered list of chains
    $self->{'chain_list'} = [];

    #read first block of alignments

    #(## ALIGNMENTS    1 -   70) line
    $line  = $text->next_line;
    ($lo, $hi) = $line =~ /^\#\#\s+ALIGNMENTS\s+(\d+)\s*-\s*(\d+)/;

    #( SeqNo  PDBNo AA STRUCTURE BP1 BP2  ACC NOCC  VAR  ....:....1..) line
    $line  = $text->next_line;
    $cut   = index($line, '..');

    ($chain1, $chaincount) = ('', 1);
     
    while (defined ($line = $text->next_line(1))) {

	if ($cut < length $line) {
	    $part2 = substr($line, $cut);
	} else {
	    #warn length($line), "($cut) $line\n";
	    $part2 = ' ' x ($hi-$lo+1);
	}
	
	$data = {};

	#chain discontinuity
	if ($line =~ /^\s*(\d+)\s+!/) {
	    (
	     $data->{'seqno'},
	     $data->{'pdbno'},
	     $data->{'aa'},
	     $data->{'structure'},
	     $data->{'bp1'},
	     $data->{'bp2'},
	     $data->{'acc'},
	     $data->{'nocc'},
	     $data->{'var'},
	    ) = ($1, '', '!', '!', 0, 0, 0, 0, 0);

	    push @{$self->{'structure'}}, $data;

	    for ($i=0; $i < $hi-$lo+1; $i++) {
		push @{$self->{'alignment'}->[$lo+$i]}, '!';
	    }

	    next;
	}

	#process the query structural information
	if ($line =~ /^
	    \s*
	    (\d+)              #SeqNo
            \s+
	    (\d+\s.)?          #PDBNo, optional, eg. '16 A'
	    \s*
	    (.)                #AA
	    \s\s               #exactly 2 spaces
	    (.........)        #STRUCTURE (9 characters)
	    (....)             #BP1
	    (.....)            #BP2
	    \s+
	    (\d+)              #ACC
	    \s+
	    (\d+)              #NOCC
	    \s+
	    (\d+)              #VAR
	    /xo) {

	    if (defined $2) {
		$self->test_args($line, $1, $2, $3, $4, $5, $6, $7, $8, $9);
	    } else {
		$self->test_args($line, $1, $3, $4, $5, $6, $7, $8, $9);
	    }

	    (
	     $data->{'seqno'},
	     $data->{'pdbno'},
	     $data->{'aa'},
	     $data->{'structure'},
	     $data->{'bp1'},
	     $data->{'bp2'},
	     $data->{'acc'},
	     $data->{'nocc'},
	     $data->{'var'},
	    ) = ($1, "$2", $3, $4, $5, $6, $7, $8, $9);


	    if ($data->{'pdbno'} =~ /^\s*$/) {
		#missing PDBNo data
		$data->{'pdbno'} = '0  ';
		$chain2 = ' ';
	    } else {
		#PDBNo is composite, eg., '16 A'
		$data->{'pdbno'} =~ /(.)$/;
		$chain2 = $1;
	    }

	    push @{$self->{'structure'}}, $data;

	    if ($chain2 eq $chain1) {
		#extend this chain
		$self->{'chain_hash'}->{$chain2}->[1] = $data->{'seqno'};
	    } else {
		#begin new chain
		#warn "|$chain1,$chain2|\n";
		$self->{'chain_hash'}->{$chain2}->[0] = $data->{'seqno'};
		$self->{'chain_hash'}->{$chain2}->[1] = $data->{'seqno'};
		push @{$self->{'chain_list'}}, $chain2;
		$chain1 = $chain2;
		$chaincount++;
	    }

# no point testing this, as maxhom miscalculates $hi sometimes.
#
#	    #process the alignments in $part2
#	    if (length $part2 > $hi-$lo+1) {
#		$self->warn("alignment columns exceed expected range: @{[length($part2)]} cf. @{[$hi-$lo+1]}");
#	    }
	    $part2 .= ' ' x ($hi - $lo + 1 - length($part2));
	    $data = [ split '', $part2 ];
	    
	    for ($i=0; $i < $hi-$lo+1; $i++) {
		push @{$self->{'alignment'}->[$lo+$i]}, $data->[$i];
	    }
	    
	    next;
	}

	last    if $line =~ /$HSSP_ALIGNMENT/o;
	last    if $line =~ /$HSSP_ALIGNMENTend/o;

	#default
	$self->warn("unknown field: $line");

    }

    return $self    unless $line =~ /$HSSP_ALIGNMENT/o;


    #read remaining blocks of alignments

    #(## ALIGNMENTS   71 -  140) line, etc
    ($lo, $hi) = $line =~ /^\#\#\s+ALIGNMENTS\s+(\d+)\s*-\s*(\d+)/;

    #( SeqNo  PDBNo AA STRUCTURE BP1 BP2  ACC NOCC  VAR  ....:....1..) line
    $line  = $text->next_line;
    $cut   = index($line, '..');

    while (defined ($line = $text->next_line(1))) {

	if ($cut < length $line) {
	    $part2 = substr($line, $cut);
	} else {
	    #warn length($line), "($cut) $line\n";
	    $part2 = ' ' x ($hi-$lo+1);
	}
	
	$data = {};

	if ($line =~ /^
	    \s*
	    (\d+)              #SeqNo
	    /xo) {

# no point testing this, as maxhom miscalculates $hi sometimes.
#
#	    #process the alignments in $part2
#	    if (length $part2 > $hi-$lo+1) {
#		$self->warn("alignment columns exceed expected range: @{[length($part2)]} cf. @{[$hi-$lo+1]}");
#	    }
	    $part2 .= ' ' x ($hi - $lo + 1 - length($part2));
	    $data = [ split '', $part2 ];
	    
	    for ($i=0; $i < $hi-$lo+1; $i++) {
		push @{$self->{'alignment'}->[$lo+$i]}, $data->[$i];
	    }
	    
	    next;
	}

	#chain discontinuity
	if ($line =~ /^\s*(\d+)\s+!/) {

	    for ($i=0; $i < $hi-$lo+1; $i++) {
		push @{$self->{'alignment'}->[$lo+$i]}, '!';
	    }
	    
	    next;
	}

	if ($line =~ /^\#\#\s+ALIGNMENTS\s+(\d+)\s*-\s*(\d+)/) {

	    ($lo, $hi) = ($1, $2); 

	    #( SeqNo  PDBNo AA STRUCTURE BP1 BP2  ACC NOCC  VAR  ....:....1..) line
	    $line  = $text->next_line;
	    $cut   = index($line, '..');

	    next;
	}

	last    if $line =~ /$HSSP_ALIGNMENTend/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'query',  $self->get_query;
    printf "$x%20s -> '%s'\n", 'dssp',   $self->get_dssp();
    #structure data [position]
    foreach my $i ($self->get_record) {
	printf "$x  %20s -> %s\n",   'seqno',      $i->{'seqno'};
	printf "$x  %20s -> '%s'\n", 'pdbno',      $i->{'pdbno'};
	printf "$x  %20s -> %s\n",   'aa',         $i->{'aa'};
	printf "$x  %20s -> '%s'\n", 'structure',  $i->{'structure'};
	printf "$x  %20s -> '%s'\n", 'bp1',        $i->{'bp1'};
	printf "$x  %20s -> '%s'\n", 'bp2',        $i->{'bp2'};
	printf "$x  %20s -> %s\n",   'acc',        $i->{'acc'};
	printf "$x  %20s -> %s\n",   'nocc',       $i->{'nocc'};
	printf "$x  %20s -> %s\n",   'var',        $i->{'var'};
    }
    #alignments[chain][rank]
    foreach my $j ($self->get_chains) {
	printf "$x%20s -> '%s'\n", "chain[$j]",  $self->get_query($j);
	for (my $i=1; $i<@{$self->{'alignment'}}; $i++) {
	    printf "$x%20s -> '%s'\n", "alignment[$j][$i]", $self->get_sequence($i, $j);
	}
    }
}

#returns list of chain names, in order of appearance in HSSP alignment
sub get_chains {
    @{$_[0]->{'chain_list'}};
}

#returns an ordered array of structure records, or just those requested in the
#order requested. record numbers are counted from 1.
sub get_record {
    my $self = shift;
    my $i;
    my @tmp = ();
    if (@_) {
	foreach $i (@_) {
	    if (defined $self->{'structure'}->[$i]) {
		push @tmp, $self->{'structure'}->[$i];
	    }
	}
    } else {
	for ($i=0; $i < @{$self->{'structure'}}; $i++) {
	    push @tmp, $self->{'structure'}->[$i];
	}
    }
    @tmp;
}

#returns the query sequence string. if a chain name is given, returns only
#that chain.
sub get_query {
    my $self = shift;
    my ($s, $i, $lo, $hi) = ('');

    #chainname supplied?
    if (@_) {
	if (exists $self->{'chain_hash'}->{$_[0]}) {
	    $lo = $self->{'chain_hash'}->{$_[0]}->[0]-1;
	    $hi = $self->{'chain_hash'}->{$_[0]}->[1];
	} else {
	    $self->warn("unknown chain: '$_[0]'");
	    return '';
	}
    } else {
	$lo = 0;
	$hi = scalar @{$self->{'structure'}};
    }

    for ($i=$lo; $i < $hi; $i++) {
	$s .= $self->{'structure'}->[$i]->{'aa'};
    }
    $s;
}

#returns the query sequence DSSP assignment string. if a chain name is
#given, returns only that chain.
sub get_dssp {
    my $self = shift;
    my ($s, $i, $lo, $hi) = ('');

    #chainname supplied?
    if (@_) {
	if (exists $self->{'chain_hash'}->{$_[0]}) {
	    $lo = $self->{'chain_hash'}->{$_[0]}->[0]-1;
	    $hi = $self->{'chain_hash'}->{$_[0]}->[1];
	} else {
	    $self->warn("unknown chain: '$_[0]'");
	    return '';
	}
    } else {
	$lo = 0;
	$hi = scalar @{$self->{'structure'}};
    }

    for ($i=$lo; $i < $hi; $i++) {
	$s .= substr($self->{'structure'}->[$i]->{'structure'}, 0, 1);
    }
    $s;
}

#returns the aligned sequence string for alignment N counting from 1. if a
#chain name is given, returns only that chain.
sub get_sequence {
    my ($self, $seqno) = (shift, shift);
    my ($s, $i, $lo, $hi) = ('');

    #chainname supplied?
    if (@_) {
	if (exists $self->{'chain_hash'}->{$_[0]}) {
	    $lo = $self->{'chain_hash'}->{$_[0]}->[0]-1;
	    $hi = $self->{'chain_hash'}->{$_[0]}->[1];
	} else {
	    $self->warn("unknown chain: '$_[0]'");
	    return '';
	}
    } else {
	$lo = 0;
	$hi = scalar @{$self->{'structure'}};
    }

    for ($i=$lo; $i < $hi; $i++) {
	my $c = $self->{'alignment'}->[$seqno]->[$i];
	if ($c eq '!') { #discontinuity
	    if ($i > 0 and $self->{'alignment'}->[$seqno]->[$i-1] eq ' ') {
		$c = ' ';
	    } elsif ($i < $hi-1 and
		     $self->{'alignment'}->[$seqno]->[$i+1] eq ' ') {
		$c = ' ';
	    }
	}
	$s .= $c;
    }
    $s;
}


###########################################################################
package NPB::Parse::Format::HSSP::PROFILE;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new { NPB::Message::warn(shift, "new() not implemented") }


###########################################################################
package NPB::Parse::Format::HSSP::INSERTION;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new { NPB::Message::warn(shift, "new() not implemented") }


###########################################################################
1;
