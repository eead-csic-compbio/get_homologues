# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: FT.pm,v 1.36 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::FT;

use NPB::Parse::Record;
use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#Exported constants
#indents for nested records
$NPB::Parse::Format::FT::indent  = 16;

#skip these many chars at front of each raw record
$NPB::Parse::Format::FT::skip    =  5;
$NPB::Parse::Format::FT::Format::FT_skip = 21;

#Generic EMBL/Genban accession number pattern
$NPB::Parse::Format::FT::ac_accession = '[A-Z][0-9]{5}';   #### CHANGE?????? 


my $MARK                  = '*';
my $Debug                 = 0;

#These recognised features were taken from EMBL release 46
my %Valid_Features = 
(
 '-'                     => 1,
 '-10_signal'            => 1,
 '-35_signal'            => 1,
 '3\'UTR'                => 1,
 '3\'clip'               => 1,
 '5\'UTR'                => 1,
 '5\'clip'               => 1,
 'CAAT_signal'           => 1,
 'CDS'                   => 1,
 'C_region'              => 1,
 'D-loop'                => 1,
 'D_segment'             => 1,
 'GC_signal'             => 1,
 'J_segment'             => 1,
 'LTR'                   => 1,
 'N_region'              => 1,
 'RBS'                   => 1,
 'STS'                   => 1,
 'S_region'              => 1,
 'TATA_signal'           => 1,
 'V_region'              => 1,
 'V_segment'             => 1,
 'allele'                => 1,
 'attenuator'            => 1,
 'conflict'              => 1,
 'enhancer'              => 1,
 'exon'                  => 1,
 'iDNA'                  => 1,
 'intron'                => 1,
 'mRNA'                  => 1,
 'mat_peptide'           => 1,
 'misc_RNA'              => 1,
 'misc_binding'          => 1,
 'misc_difference'       => 1,
 'misc_feature'          => 1,
 'misc_recomb'           => 1,
 'misc_signal'           => 1,
 'misc_structure'        => 1,
 'modified_base'         => 1,
 'mutation'              => 1,
 'old_sequence'          => 1,
 'polyA_signal'          => 1,
 'polyA_site'            => 1,
 'precursor_RNA'         => 1,
 'prim_transcript'       => 1,
 'primer_bind'           => 1,
 'promoter'              => 1,
 'protein_bind'          => 1,
 'rRNA'                  => 1,
 'rep_origin'            => 1,
 'repeat_region'         => 1,
 'repeat_unit'           => 1,
 'satellite'             => 1,
 'scRNA'                 => 1,
 'sig_peptide'           => 1,
 'snRNA'                 => 1,
 'source'                => 1,
 'stem_loop'             => 1,
 'tRNA'                  => 1,
 'terminator'            => 1,
 'transit_peptide'       => 1,
 'unsure'                => 1,
 'variation'             => 1,
 'virion'                => 1,
);


sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);
    
    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::FT::skip);

    while (defined ($line = $text->next_line)) {

	if ($line =~ /^(\S+)/) {
	    if (exists $Valid_Features{$1}) {
		$text->scan_nest($NPB::Parse::Format::FT::indent, $1);
		next;
	    }
	}
	
	#default
	$line = substr($line, 0, $NPB::Parse::Format::FT::indent);
	$line =~ tr/ //d;
	$text->scan_nest($NPB::Parse::Format::FT::indent, "$line");
        $self->warn("unknown field: $line");
    }
    $self;
}

#Override NPB::Parse::Record::get_class() which returns the Record hierarchy
#of $self (which could involve Embl or Genbank, etc), to return the generic
#feature table class defined here. Used by NPB::Parse::Record::add_object().
sub get_class {
    'NPB::Parse::Format::FT';
}


###########################################################################
package NPB::Parse::Feature;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    shift;               #discard system supplied type
    my $type = shift;    #use caller supplied type
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);
    
    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::FT::Format::FT_skip);

    $self->{'location'}     = '';
    $self->{'qualifiers'}   = '';
    $self->{'p_location'}   = [];
    $self->{'p_qualifiers'} = {};

    #location lines: scan until first qualifier
    $record = $text->scan_until('^\/');
    $self->parse_location($record);

    #qualifier lines: scan until end-of-record
    $record = $text->scan_lines(0);
    $self->parse_qualifiers($record);

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",   'location',   $self->{'location'};
    printf "$x%20s -> %s\n",   'qualifiers', $self->{'qualifiers'};
    printf "$x%20s -> [%s]\n", 'p_location', join(",",@{$self->{'p_location'}});
    printf "$x%20s:\n", 'p_qualifiers';
    foreach my $i (sort keys %{$self->{'p_qualifiers'}}) {
	printf("$x$x%20s -> [%s]\n",
               $i, join(',', @{$self->{'p_qualifiers'}->{$i}}));
    }
}

sub parse_location {
    my ($self, $loc) = @_;
    my @parse;
    my $make;

    #strip newlines and blanks
    $loc =~ tr/ \t\n//d;

    $self->{'location'} = $loc;

    if ((@parse = $self->_parse_location($loc)) < 1) {
	$self->warn("unrecognised operator in: $loc");
	return;
    }
    $make = $self->_make_location(@parse);

    warn "PARSE= [", join(', ', @parse), "]\n"    if $Debug;
    warn "ORIG=  [$loc]\n"                        if $Debug;
    warn "MADE=  [$make]\n"                       if $Debug;

    if ($loc ne $make) {
	$self->warn("parse_location() original and reconstructed parse trees differ:\n    orig=  [$loc]\n    parse= [$make]");
	return;
    }

    $self->{'p_location'} = \@parse;
}

sub parse_qualifiers {
    my ($self, $text) = @_;
    my ($qualifier, $value);

    $self->{'qualifiers'} = $text;

    while ($text) {

	#read qualifier name from single line
	if ($text =~ /^(\/[^\/= ]+)/) {
	    ($qualifier, $text) = ($1, $');
	} else {
	    $self->warn("unmatched qualifier: $text");
	    last;
	}
	
	#read optional value, possibly over multiple lines, stripping
	#newlines where present.

	#quoted string
	if ($text =~ /^\s*=\s*(\".*)/s) {
	    #print STDERR "LOOKING AT: ($qualifier) quoted string\n";
	    ($value, $text) = _match_value($1);
	    
	    if      ($qualifier eq '/note') {
		$value = NPB::Parse::Record::strip_english_newlines($value);
	    } elsif ($qualifier eq '/function') {
		$value = NPB::Parse::Record::strip_english_newlines($value);
	    } elsif ($qualifier eq '/product') {
		$value = NPB::Parse::Record::strip_english_newlines($value);
	    } else {
		$value =~ tr/\n//d;
		$value =~ s/^\"//;
		$value =~ s/\"$//;
	    }

	    #print STDERR "READ VALUE: $value\n";
	    push @{$self->{'p_qualifiers'}->{$qualifier}}, $value;
	    next;
	}

	#bracketed list
	if ($text =~ /^\s*=\s*(\(.*)/s) {
	    #print STDERR "LOOKING AT: ($qualifier) bracketed list\n";
	    ($value, $text) = _match_value($1);
	    $value =~ tr/\n//d;
	    #print STDERR "READ VALUE: $value\n";
	    push @{$self->{'p_qualifiers'}->{$qualifier}}, $value;
	    next;
	}

	#single value
	if ($text =~ /^\s*=\s*(.*)/s) {
	    #print STDERR "LOOKING AT: ($qualifier) single value\n";
	    ($value, $text) = _match_value($1);
	    $value =~ tr/\n//d;
	    #print STDERR "READ VALUE: $value\n";
	    push @{$self->{'p_qualifiers'}->{$qualifier}}, $value;
	    next;
	}

	#simple flag
	if ($text =~ /^\s*(.*)/s) {
	    #print STDERR "LOOKING AT: ($qualifier) simple flag\n";
	    ($value, $text) = _match_value($1);
	    #print STDERR "READ VALUE: $value\n";
	    push @{$self->{'p_qualifiers'}->{$qualifier}}, $value;
	    next;
	}
	
	#default
	$self->warn("unmatched value: $text");
    }
}

sub _match_value {
    my $text = shift;

    #test paired " ignoring double "" and contained newlines
    if ($text =~ /^
	\n*
	(
	 \"
	 \n*
	 ([^\"] | \"\n*\")*
	 \n*
	 \"
	 \n*
	)/sxo) {
	return ($1, $');
    } 

    #test paired non-nested ( ), ignoring contained newlines
    if ($text =~ /^
	\n*
	(
	 \(
	 \n*
	 [^\)]*
	 \n*
	 \)
	 \n*
	)/sxo) {
	return ($1, $');
    }

    #test single value, or null value
    if ($text =~ /^
	\n*
	(
	 [^\/]*
	 \n*
	)/sxo) {
	return ($1, $');
    }

    ('', $text)
}	

sub _parse_location {
    my ($self, $loc) = @_;
    my @parse = ();
    my $op;

    warn "_parse_location: $loc\n"    if $Debug;

    #MULTIPLE OPERATORS '(join|order|group|one-of)'
    if ($loc =~ /^(join|order|group|one-of)\((.*)\)$/o) {
	$op = $1;
	if ((@parse = $self->_parse_list($2, ',')) > 0) {
	    return ($MARK, @parse, "M($op)");
	}
	return ();
    }

    #BINARY 'replace'
    if ($loc =~ /^replace\((.*),(.*)\)$/o) {
	if ((@parse = $self->_parse_binary($1, $2)) > 0) {
	    return ($MARK, @parse, 'B(replace)');
	}
	return ();
    }

    #UNARY 'complement'
    if ($loc =~ /^complement\((.*)\)$/o) {
	if ((@parse = $self->_parse_location($1)) > 0) {
	    return ($MARK, @parse, 'U(complement)');
	}
	return ();
    }

    #string constant
    if ($loc =~ /^"([^\"]*)"$/o) {
	return ($MARK, "K(\"$1\")");
    }

    #label
    if ($loc =~ /^([-\'\w]+)$/o) {
	return ($MARK, "L($1)");
    }

    #cross-reference
    if ($loc =~ /^($NPB::Parse::Format::FT::ac_accession):(.*)$/o) {
	$op = $1;
	if ((@parse = $self->_parse_location($2)) > 0) {
	    return ($MARK, @parse, "X($op)");
	}
	return ();
    }

    #BINARY OPERATORS '(..|.|^)'
    #  (A.B)..(C.D)
    #  (A.B)..C
    #  (A..(B.C)
    #  (A..B)
    #  A..B

    #left and right bracketted expression around operator: (900.901)..(906.907)
    if ($loc =~ /^(\(.*\))(\.\.|\.|\^)(\(.*\))$/o) {
	$op = $2;
	if ((@parse = $self->_parse_binary($1, $3)) > 0) {
	    return ($MARK, @parse, "I($op)");
	}
	return ();
    }

    #left bracketted expression before operator: (900.901)..906
    if ($loc =~ /^(\(.*\))(\.\.|\.|\^)([^(].*[^)])$/o) {
	$op = $2;
	if ((@parse = $self->_parse_binary($1, $3)) > 0) {
	    return ($MARK, @parse, "I($op)");
	}
	return ();
    }

    #right bracketted expression after operator: 900..(901.906)
    if ($loc =~ /^([^(].*[^).^])(\.\.|\.|\^)(\(.*\))$/o) {
	$op = $2;
	if ((@parse = $self->_parse_binary($1, $3)) > 0) {
	    return ($MARK, @parse, "I($op)");
	}
	return ();
    }

    #bracketted simple operator: (901..906)
    if ($loc =~ /^\((.*[^.^])(\.\.|\.|\^)(.*)\)$/o) {
	$op = $2;
	if ((@parse = $self->_parse_binary($1, $3)) > 0) {
	    return ($MARK, @parse, "P($op)");
	}
	return ();
    }

    #simple operator: 901..906
    if ($loc =~ /^(.*[^.^])(\.\.|\.|\^)(.*)$/o) {
	$op = $2;
	if ((@parse = $self->_parse_binary($1, $3)) > 0) {
	    return ($MARK, @parse, "I($op)");
	}
	return ();
    }

    #failed
    $self->warn("unrecognised operator in: $loc");

    return ();
}

sub _parse_list {
    my ($self, $loc, $pat) = @_;
    my (@parse, $expr, @result) = (());

    warn "_parse_list:     $loc\n"    if $Debug;

    foreach $expr ($self->_expand_one_level($loc, $pat)) {
	if ((@result = $self->_parse_location($expr)) < 1) {
	    return ();
	}
	push @parse, @result;
    }
    return @parse;
}

sub _parse_binary {
    my ($self, $loc1, $loc2) = @_;
    my (@parse, $expr, @result) = (());

    warn "_parse_binary:   $loc1 / $loc2\n"    if $Debug;

    foreach $expr ($loc1, $loc2) {

	#bracketted?
	if ($expr =~ /^(\(.*\))$/o) {
	    if ((@result = $self->_parse_location($1)) < 1) {
		return ();
	    }
	    push @parse, @result;
	    next;
	}

	#numeric item with < or > qualifier
	if ($expr =~ /^([<>])(\d+)$/o) {
	    push @parse, ($MARK, $2, "S($1)");
	    next;
	}
	
	#numeric item
	if ($expr =~ /^(\d+)$/o) {
	    push @parse, $1;
	    next;
	}

	#string constant (used by 'replace' operator)
	if ($expr =~ /^("[^\"]*")/o) {
	    if ((@result = $self->_parse_location($1)) < 1) {
		return ();
	    }
	    push @parse, @result;
	    next;
	}

	#failed
	$self->warn("unrecognised operator in: $expr");
	return ();
    }
    return @parse;
}

#expand a string composed of a $pattern separated list by one level, leaving 
#nested lists within parenthetic pairs. return the list so formed.
sub _expand_one_level {
    my ($self, $string, $pattern) = @_;
    my ($depth, $c, $last, @list) = (0, 0, 0, ());
    my ($plen, $slen) = (length($pattern), length($string));

    warn("_expand_one_level: entry= $string\n")    if $Debug;

    for ($c=0; $c < $slen; $c++) {
	if (substr($string, $c, 1) eq "(") {
	    $depth++;
	} elsif (substr($string, $c, 1) eq ")") {
	    $depth--;
	} elsif (substr($string, $c, $plen) eq $pattern) {
	    if ($depth == 0) {
		push @list, substr($string, $last, $c-$last);
		$last = $c + $plen;
		warn("_expand_one_level: push=  ", join("/", @list), "\n")    if $Debug;
	    }
	    $c += $plen - 1;
	}
    }
    if ($depth == 0 and $last < $slen) {
	push @list, substr($string, $last, $c-$last);
    } elsif ($depth) {
	$self->warn("unmatched parenthesis in: $string");
    }
    warn("_expand_one_level: exit=  ", join("/", @list), "\n")    if $Debug;
    return @list;
}

sub _make_location {
    my ($self, @tree) = @_;
    my ($token, @tmp);
    my @opstack  = ();
    my @argstack = ();

    while ($token = pop @tree) {
	if ($token =~ /^[SUBPIMKLX]/) {
	    #operator
	    push @opstack, $token;
	} elsif ($token ne $MARK) {
	    #argument
	    push @argstack, $token;
	} else {
	    #mark
	    $token = pop @opstack;
	    if      ($token =~ /^S\((.*)\)/) {
		$token = $1;
		@tmp = (pop @argstack);
		push @argstack, $token . join('', @tmp);
		warn "SIGN:   $token, @tmp  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^U\((.*)\)/) {
		$token = $1;
		@tmp = (pop @argstack);
		push @argstack, $token . '(' . join('', @tmp) . ')';
		warn "UNARY:  $token, @tmp  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^B\((.*)\)/) {
		$token = $1;
		@tmp = (pop @argstack); push @tmp, pop @argstack;
		push @argstack, $token . '(' . join(',', @tmp) . ')';
		warn "BINARY: $token, @tmp  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^I\((.*)\)/) {
		$token = $1;
		@tmp = (pop @argstack); push @tmp, pop @argstack;
		push @argstack, $tmp[0] . $token . $tmp[1];
		warn "INFIX:  $token, @tmp  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^P\((.*)\)/) {
		$token = $1;
		@tmp = (pop @argstack); push @tmp, pop @argstack;
		push @argstack, '(' . $tmp[0] . $token . $tmp[1] . ')';
		warn "PAREN:  @tmp  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^M\((.*)\)/) {
		push @opstack, $1;
		@tmp = (); while ($token = pop @argstack) {
		    if ($token eq $MARK) {
			push @argstack, $token;
			last;
		    }
		    push @tmp, $token;
		}
		$token = pop @opstack;
		push @argstack, $token . '(' . join(',', @tmp) . ')';
		warn "MULTI:  $token, @tmp  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^K\((.*)\)/) {
		$token = $1;
		push @argstack, $token;
		warn "STRING: $token  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^L\((.*)\)/) {
		$token = $1;
		push @argstack, $token;
		warn "LABEL:  $token  =>  @argstack\n"    if $Debug;
	    } elsif ($token =~ /^X\((.*)\)/) {
		push @opstack, $1;
		@tmp = (); while ($token = pop @argstack) {
		    if ($token eq $MARK) {
			push @argstack, $token;
			last;
		    }
		    push @tmp, $token;
		}
		$token = pop @opstack;
		push @argstack, $token . ':' . join(',', @tmp);
		warn "XREF:  $token, @tmp  =>  @argstack\n"    if $Debug;
	    }
	}
    }
    if (@opstack or @argstack != 1) {
	$self->warn("bad stacks: op=[", join(",", @opstack), "] arg=[", join(",", @argstack), "] tree=[", join(",", @tree), "]");
    }
    join('', @argstack);
}


###########################################################################
package NPB::Parse::Format::FT::CDS;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Feature);

sub new {
    new NPB::Parse::Feature(@_);
}


###########################################################################
package NPB::Parse::Format::FT::tRNA;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Feature);

sub new {
    new NPB::Parse::Feature(@_);
}


###########################################################################
package NPB::Parse::Format::FT::rRNA;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Feature);

sub new {
    new NPB::Parse::Feature(@_);
}


###########################################################################
package NPB::Parse::Format::FT::exon;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Feature);

sub new {
    new NPB::Parse::Feature(@_);
}


###########################################################################
1;
