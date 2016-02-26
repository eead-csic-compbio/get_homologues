# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Sequence.pm,v 1.15 2015/06/14 17:09:03 npb Exp $

###########################################################################
package Bio::MView::Sequence;
use strict;

use vars qw($WARNCLASH $OVERWRITE $HARDFASTA);

my $Find_Pad = '[-._~]';  #input terminal gap characters
my $Find_Gap = '[-._~]';  #input internal gap characters
my $Find_Spc = '\s';      #input whitespace character
my $Find_Fs1 = '\/';      #input frameshift mark '/' from fast[xy], tfast[axy]
my $Find_Fs2 = '\\\\';    #input frameshift mark '\' from fast[xy], tfast[axy]

my $Text_Pad = '.';       #default output terminal gap
my $Text_Gap = '-';       #default output internal gap
my $Text_Spc = ' ';       #default output whitespace
my $Text_Fs1 = '/';       #default output frameshift1
my $Text_Fs2 = '\\';      #default output frameshift2

my $Mark_Pad = "\001";    #encoded terminal gap
my $Mark_Gap = "\002";    #encoded internal gap
my $Mark_Spc = "\003";    #encoded whitespace
my $Mark_Fs1 = "\004";    #encoded frameshift
my $Mark_Fs2 = "\005";    #encoded frameshift

$WARNCLASH = 0;    #set to 1 to enable reporting of clashing symbols
$OVERWRITE = 0;    #set to 1 to enable clash overwrite by newest symbol
$HARDFASTA = 0;    #relax length rules for UV. FASTA input

#All range numbers count from 1, set to 0 when undefined.
#
#Externally, columns are numbered (increasing) left to right from 1.
#
#Internally, each row is numbered in relation to the reference sequence
#against which all other rows are being assembled.
#
#'lo' and 'hi' for any row give the lowest and highest (by magnitude) sequence
#positions in the reference sequence for that row including terminal gaps.
#
#'reflo' and 'refhi' for any row give the starting and stopping positions in
#the reference sequence in the orientation of the reference seqence for that
#row excluding terminal gaps.
#
#'pfxlen' and 'sfxlen' give the lengths of the terminal gaps for that row:
#
# Forwards query:

# column  12345
# +query  23456 (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 2,6, 0,0)
# +hit    -23-- (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 3,4, 1,2)
#
# column  12345
# +query  23456 (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 2,6, 0,0)
# -hit    -32-- (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 3,4, 1,2)
#
# Reverse query:
#
# column  12345
# -query  65432 (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 6,2, 0,0)
# +hit    -23-- (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 5,4, 1,2)
#
# column  12345
# -query  65432 (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 6,2, 0,0)
# -hit    -32-- (lo,hi,	reflo,refhi, pfxlen,sfxlen) = (2,6, 5,4, 1,2)
#
sub new {
    my $type = shift;
    my $self = {};

    $self->{'seq'}    	= {};   	#sequence hash of char

    $self->{'lo'}     	= 0;    	#absolute start of string
    $self->{'hi'}     	= 0;    	#absolute end of string

    $self->{'pfxlen'} 	= 0;    	#length of non-sequence prefix
    $self->{'sfxlen'} 	= 0;    	#length of non-sequence suffix

    $self->{'reflo'} 	= 0;    	#first position of sequence data
    $self->{'refhi'} 	= 0;    	#last  position of sequence data

    $self->{'f1'}  	= undef;    	#from label of sequence
    $self->{'t1'}   	= undef;    	#to   label of sequence

    $self->{'f2'}  	= undef;    	#2nd from label of sequence (optional)
    $self->{'t2'}   	= undef;    	#2nd to   label of sequence (optional)

    $self->{'fs'}       = 0;            #has frameshifts

    $self->{'find_pad'} = $Find_Pad;
    $self->{'find_gap'} = $Find_Gap;
    $self->{'find_spc'} = $Find_Spc;
    $self->{'find_fs1'} = $Find_Fs1;
    $self->{'find_fs2'} = $Find_Fs2;

    $self->{'text_pad'} = $Text_Pad;
    $self->{'text_gap'} = $Text_Gap;
    $self->{'text_spc'} = $Text_Spc;
    $self->{'text_fs1'} = $Text_Fs1;
    $self->{'text_fs2'} = $Text_Fs2;

    bless $self, $type;

    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

sub print {
    sub _format {
	my ($self, $k, $v) = @_;
	$v = 'undef' unless defined $v;
	$v = "$v [@{[$self->string]}]" if $k eq 'seq';
	$v = "'$v'" if $v =~ /^\s*$/;
	return sprintf("  %-15s => %s\n", $k, $v)
    }
    my $self = shift;
    warn "$self\n";
    map { warn $self->_format($_, $self->{$_}) } sort keys %$self;
    $self;
}

sub encode {
    my ($self, $s) = @_;

    #leading non-sequence characters
    while ($$s =~ s/^($Mark_Pad*)$self->{'find_gap'}/$1$Mark_Pad/g) {}

    #trailing non-sequence characters
    while ($$s =~ s/$self->{'find_gap'}($Mark_Pad*)$/$Mark_Pad$1/g) {}

    #internal gap characters
    $$s =~ s/$self->{'find_gap'}/$Mark_Gap/g;

    #internal spaces
    $$s =~ s/$self->{'find_spc'}/$Mark_Spc/g;

    #frameshift '/'
    $$s =~ s/$self->{'find_fs1'}/$Mark_Fs1/g;

    #frameshift '\'
    $$s =~ s/$self->{'find_fs2'}/$Mark_Fs2/g;

    $self;
}

sub string {
    my $s = $_[0]->_substr($_[0]->{'lo'}, $_[0]->{'hi'});
    $s =~ s/$Mark_Pad/$_[0]->{'text_pad'}/g;
    $s =~ s/$Mark_Gap/$_[0]->{'text_gap'}/g;
    $s =~ s/$Mark_Spc/$_[0]->{'text_spc'}/g;
    $s =~ s/$Mark_Fs1/$_[0]->{'text_fs1'}/g;
    $s =~ s/$Mark_Fs2/$_[0]->{'text_fs2'}/g;
    $s;
}

sub sequence {
    my $s = $_[0]->_substr($_[0]->{'lo'}, $_[0]->{'hi'});
    $s =~ s/$Mark_Pad//og;
    $s =~ s/$Mark_Gap//og;
    $s =~ s/$Mark_Spc//og;
    $s =~ s/$Mark_Fs1//og;
    $s =~ s/$Mark_Fs2//og;
    $s;
}

# sub length {
#     my $self = shift;
#     my $claim = $self->lablen;
#     my $true = $self->trulen;
#     warn "$self: length: true: $true  claim: $claim\n";
#     warn "T= [@{[$self->string]}]\n";
#     $claim;
# }

#return actual string length; includes special characters
sub length { $_[0]->trulen }
#sub length { $_[0]->lablen }

#return pure sequence string length; excludes special characters
sub seqlen { CORE::length $_[0]->sequence }

#return actual string length ; includes special characters
sub trulen { CORE::length $_[0]->string }

#return string length claimed by sequence labels; includes special characters
sub lablen {
    return $_[0]->{'hi'} - $_[0]->{'lo'} + 1  if $_[0]->{'lo'} > 0;
    return 0;
}

#return real oriented begin/end positions
sub lo { $_[0]->{'lo'} }
sub hi { $_[0]->{'hi'} }

#return original from/to labels (default, or eg., search query)
sub fromlabel1 { $_[0]->{'f1'} }
sub tolabel1   { $_[0]->{'t1'} }

#return original from/to labels, second set (eg., search hit)
sub fromlabel2 { $_[0]->{'f2'} }
sub tolabel2   { $_[0]->{'t2'} }

#return internal start/stop position of sequence wrt reference
sub reflo { $_[0]->{'reflo'} }
sub refhi { $_[0]->{'refhi'} }

#return lengths of leading/trailing terminal gap regions
#sub leader  { my $n = $_[0]->{'reflo'} - $_[0]->{'lo'}; $n > -1 ? $n : 0 }
#sub trailer { my $n = $_[0]->{'hi'} - $_[0]->{'refhi'}; $n > -1 ? $n : 0 }
sub leader  { $_[0]->{'pfxlen'} }
sub trailer { $_[0]->{'sfxlen'} }

#return count of frameshifts
sub frameshifts { $_[0]->{'fs'} }

sub set_find_pad { $_[0]->{'find_pad'} = $_[1] }
sub set_find_gap { $_[0]->{'find_gap'} = $_[1] }
sub set_find_spc { $_[0]->{'find_spc'} = $_[1] }
sub set_find_fs1 { $_[0]->{'find_fs1'} = $_[1] }
sub set_find_fs2 { $_[0]->{'find_fs2'} = $_[1] }

sub set_pad { $_[0]->{'text_pad'} = $_[1] }
sub set_gap { $_[0]->{'text_gap'} = $_[1] }
sub set_spc { $_[0]->{'text_spc'} = $_[1] }
sub set_fs1 { $_[0]->{'text_fs1'} = $_[1] }
sub set_fs2 { $_[0]->{'text_fs2'} = $_[1] }

sub get_pad { $_[0]->{'text_pad'} }
sub get_gap { $_[0]->{'text_gap'} }
sub get_spc { $_[0]->{'text_spc'} }
sub get_fs1 { $_[0]->{'text_fs1'} }
sub get_fs2 { $_[0]->{'text_fs2'} }

sub set_range {
    my ($self, $lo, $hi) = @_;
    die "$self: range values in wrong order ($lo, $hi)\n"    if $lo > $hi;
    ($self->{'lo'}, $self->{'hi'}) = ($lo, $hi);
    $self;
}

sub is_reversed {0}

sub reverse {
    no strict qw(subs);
    bless $_[0], Bio::MView::Reverse_Sequence;
}

#Horribly inefficient search for all matches of a regexp in the possibly
#gapped sequence, ie., regexp matches can span gapchars. Returns a reference
#to a list of matched non-gap sequence positions (indexing from 1).
sub findall {
    my ($self, $blocks, $mapsize) = @_;
    my $list = [];
    my $sequence = $self->sequence;
    my $seqlen = CORE::length $sequence;
    my @gapped = split(//, $self->string);
    my $gapchar = $self->{text_gap};

    #warn "findall(): @{[$self->string]}\n";
    #warn "findall(): BLOCK=$pattern\n";

    my $end = scalar @$blocks;

    for (my $blocknum=0; $blocknum<$end; $blocknum++) {
	#warn "findall(): $blocknum, $blocks->[$blocknum]\n";
	my $blockid = chr(ord('A') + $blocknum % $mapsize);
	my $pattern = $blocks->[$blocknum];
	study $pattern;

	for (my ($i,$j) = (0,0); $j<$seqlen; $i++, $j++) {
	    while ($gapped[$i] eq $gapchar) {
		$i++;  #sync gapped/ungapped string indices
	    }
	    my $s = substr($sequence, $j);
	    if ($s =~ /^($pattern)/i) {  #case insensitive
		my @match = split(//, $1);
		#warn "\nfindall(): $i,$j [@match] $s\n";
		for (my ($k,$l) = ($i,0); $l<@match; $k++, $l++) {
		    while ($gapped[$k] eq $gapchar) {
			$k++;  #absorb gaps
		    }
		    push @$list, [$blockid, $k+1];
		}
	    }
	}
    }
    return $list;
}

#input each argument frag as [string, from, to], where from <= to
sub insert {
    my $self = shift;
    my ($string, $frag, $len, $i, $p, $c, $state);

    $state = 0;

    foreach $frag (@_) {
	die "insert(+) wrong direction ($frag->[1], $frag->[2])\n"
	    if $frag->[1] > $frag->[2];

	$string = ${ $frag->[0] };
        #warn "frg(+)=@$frag\n";
        #warn "frg(+)=$frag->[1], $frag->[2] [$string]\n";

	$self->encode(\$string);
	#warn "frg(+)=$frag->[1], $frag->[2] [$string]\n";

	$len = CORE::length $string;

        #update sequence ranges
        my $lo = $frag->[1];
	my $hi = $frag->[2];

        $self->{'lo'} = $lo  if $lo < $self->{'lo'} or $self->{'lo'} == 0;
        $self->{'hi'} = $hi  if $hi > $self->{'hi'} or $self->{'hi'} == 0;

        #update sequence labels: optional first pair
        if (defined $frag->[3] and $frag->[3] != 0) {
	    $self->{'f1'} = $frag->[3]   unless defined $self->{'f1'};
	    $self->{'t1'} = $frag->[4]   unless defined $self->{'t1'};
	    my $f = $frag->[3];
	    my $t = $frag->[4];
	    if ($f < $t or $self->{'f1'} < $self->{'t1'}) {
		#label pair orientated positive
		$self->{'f1'} = $f  if $f < $self->{'f1'};
		$self->{'t1'} = $t  if $t > $self->{'t1'};
	    } else {
		#label pair orientated negative
		$self->{'f1'} = $f  if $f > $self->{'f1'};
		$self->{'t1'} = $t  if $t < $self->{'t1'};
	    }
	}

        #update sequence labels: optional second pair
        if (defined $frag->[5] and $frag->[5] != 0) {
	    $self->{'f2'} = $frag->[5]   unless defined $self->{'f2'};
	    $self->{'t2'} = $frag->[6]   unless defined $self->{'t2'};
	    my $f = $frag->[5];
	    my $t = $frag->[6];
	    if ($f < $t or $self->{'f2'} < $self->{'t2'}) {
		#label pair orientated positive
		$self->{'f2'} = $f  if $f < $self->{'f2'};
		$self->{'t2'} = $t  if $t > $self->{'t2'};
	    } else {
		#label pair orientated negative
		$self->{'f2'} = $f  if $f > $self->{'f2'};
		$self->{'t2'} = $t  if $t < $self->{'t2'};
	    }
	}

	#populate sequence array
	for ($i=0; $i < $len; $i++) {

	    $c = substr($string, $i, 1);

	    #begin/end
	    if ($c eq $Mark_Pad) {
		if ($state == 0) {
		    $self->{'pfxlen'}++;
		} elsif ($state > 0) {
		    $self->{'sfxlen'}++;
		    $state = 2;
		}
		next;
	    }

	    #middle
	    $state = 1    if $state < 1;

	    #skip gaps
	    next    if $c eq $Mark_Gap;

	    $p = $lo + $i;

	    #warn "insert(+): $p/($self->{'lo'},$self->{'hi'}) = $i/$len\t[$c]\n";

	    #store other text, including Mark_Spc or Mark_Fs[12] symbols
	    if (exists $self->{'seq'}->{$p}) {
		next if $self->{'seq'}->{$p} eq $c;
		warn "${self}::insert(+): assembly clash: [$p] = '$$self->{'seq'}->{$p}' -> '$c'\n"
		    if $Bio::MView::Sequence::WARNCLASH;
		next unless $Bio::MView::Sequence::OVERWRITE;
		#unrecord old frameshift
		$self->{'fs'}--  if $self->{'seq'}->{$p} eq $Mark_Fs1;
		$self->{'fs'}--  if $self->{'seq'}->{$p} eq $Mark_Fs2;
	    }
	    #record frameshift
	    $self->{'fs'}++  if $c eq $Mark_Fs1;
	    $self->{'fs'}++  if $c eq $Mark_Fs2;

	    #write the character
	    $self->{'seq'}->{$p} = $c;
	}

	
	if ($Bio::MView::Sequence::HARDFASTA and
	    $len != $hi-$lo+1 and $self->{'fs'} < 1) {
	    warn "${self}::insert(+) unframeshifted length mismatch ($len != $hi-$lo+1)\n";
	}
	#warn "insert(+): $self->{'lo'} $self->{'hi'}\n";
    }

    #adjust prefix/suffix positions given new lengths
    $self->{'reflo'} = $self->{'lo'} + $self->{'pfxlen'};
    $self->{'refhi'} = $self->{'hi'} - $self->{'sfxlen'};

    #$self->print;

    $self;
}

sub _substr {
    my ($self, $start, $stop) = @_;

    return ''  if $start < 1 or $stop < 1;   #negative range args
    return ''  unless $self->{'lo'} > 0;     #empty
    return ''  if $stop  < $self->{'lo'};    #missed (too low)
    return ''  if $start > $self->{'hi'};    #missed (too high)

    $start = $self->{'lo'}  if $start < $self->{'lo'};
    $stop  = $self->{'hi'}  if $stop  > $self->{'hi'};
    #warn "_substr(+,lo,hi): $start, $stop";
    #warn "_substr(+,beg,end): $self->{'reflo'}, $self->{'refhi'}";
    $stop++;

    my $s = '';

    for (my $i = $start; $i < $stop; $i++) {
	if ($i < $self->{'reflo'} or $i > $self->{'refhi'}) {
	    $s .= $Mark_Pad;
	    next;
	}
	#warn "_substr(+): $i [$self->{'seq'}->{$i}]";
	if (exists $self->{'seq'}->{$i}) {
	    $s .= $self->{'seq'}->{$i};
	} elsif ($i > $self->{'reflo'} and $i < $self->{'refhi'}) {
	    #We never write a gap at the edge. Frameshifts upset the
	    #numbering so that [begin,end] may exceed the true range.
	    #warn "_substr(+): gap: [$i] ($self->{lo},$self->{hi}), ($self->{reflo},$self->{refhi})\n";
	    $s .= $Mark_Gap;
	}
    }
    $s;
}

sub raw {
    my ($self, $col, $p) = @_;

    $p = $self->{'lo'} + $col - 1;

    return ''  if $p < $self->{'lo'} or $p > $self->{'hi'};

    return $Mark_Pad
	if $p < $self->{'reflo'} or $p > $self->{'refhi'};
    return $self->{'seq'}->{$p}
        if exists $self->{'seq'}->{$p};
    return $Mark_Gap;
}

sub col {
    my ($self, $col, $p) = @_;

    $p = $self->{'lo'} + $col - 1;

#    warn("col(+): $col [", 
#	 exists $self->{'seq'}->{$col} ? $self->{'seq'}->{$col} : '',
#	 "]=> $p [", 
#	 exists $self->{'seq'}->{$p} ? $self->{'seq'}->{$p} : '',
#	 "]\n");

    return ''  if $p < $self->{'lo'} or $p > $self->{'hi'};

    return $self->{'text_pad'} if $p < $self->{'reflo'} or $p > $self->{'refhi'};
    return $self->{'text_gap'}  unless exists $self->{'seq'}->{$p};
    return $self->{'text_spc'} if $self->{'seq'}->{$p} eq $Mark_Spc;
    return $self->{'text_fs1'} if $self->{'seq'}->{$p} eq $Mark_Fs1;
    return $self->{'text_fs2'} if $self->{'seq'}->{$p} eq $Mark_Fs2;
    return $self->{'seq'}->{$p};
}

sub is_char {
    $_[1] ne $Mark_Pad and $_[1] ne $Mark_Gap and $_[1] ne $Mark_Spc
	and $_[1] ne $Mark_Fs1 and $_[1] ne $Mark_Fs2;
}
sub is_non_char     { ! $_[0]->is_char($_[1]) }
sub is_space        { $_[1] eq $Mark_Spc }
sub is_frameshift   { $_[1] eq $Mark_Fs1 or $_[1] eq $Mark_Fs2 }
sub is_terminal_gap { $_[1] eq $Mark_Pad }
sub is_internal_gap { $_[1] eq $Mark_Gap }
sub is_gap          { $_[1] eq $Mark_Pad or $_[1] eq $Mark_Gap }

###########################################################################
package Bio::MView::Reverse_Sequence;

use Bio::MView::Sequence;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Sequence);

sub is_reversed {1}

sub reverse {
    no strict qw(subs);
    bless $_[0], Bio::MView::Sequence;
}

#input each argument frag as [string, from, to], where from >= to
sub insert {
    my $self = shift;
    my ($string, $frag, $len, $i, $p, $c, $state);

    $state = 0;

    foreach $frag (@_) {
	die "insert(-) wrong direction ($frag->[1], $frag->[2])\n"
	    if $frag->[2] > $frag->[1];

	$string = ${ $frag->[0] };
        #warn "frg(-)=@$frag\n";
        #warn "frg(-)=$frag->[1], $frag->[2] [$string]\n";

	$self->encode(\$string);
	#warn "frg(-)=$frag->[1], $frag->[2] [$string]\n";

	$len = CORE::length $string;
	
        #update sequence ranges REVERSE
        my $lo = $frag->[2];
	my $hi = $frag->[1];

        $self->{'lo'} = $lo  if $lo < $self->{'lo'} or $self->{'lo'} == 0;
        $self->{'hi'} = $hi  if $hi > $self->{'hi'} or $self->{'hi'} == 0;

        #update sequence labels: optional first pair
        if (defined $frag->[3] and $frag->[3] != 0) {
	    $self->{'f1'} = $frag->[3]   unless defined $self->{'f1'};
	    $self->{'t1'} = $frag->[4]   unless defined $self->{'t1'};
	    my $f = $frag->[3];
	    my $t = $frag->[4];
	    if ($f < $t or $self->{'f1'} < $self->{'t1'}) {
		#label pair orientated positive
		$self->{'f1'} = $f  if $f < $self->{'f1'};
		$self->{'t1'} = $t  if $t > $self->{'t1'};
	    } else {
		#label pair orientated negative
		$self->{'f1'} = $f  if $f > $self->{'f1'};
		$self->{'t1'} = $t  if $t < $self->{'t1'};
	    }
	}

        #update sequence labels: optional second pair
        if (defined $frag->[5] and $frag->[5] != 0) {
	    $self->{'f2'} = $frag->[5]   unless defined $self->{'f2'};
	    $self->{'t2'} = $frag->[6]   unless defined $self->{'t2'};
	    my $f = $frag->[5];
	    my $t = $frag->[6];
	    if ($f < $t or $self->{'f2'} < $self->{'t2'}) {
		#label pair orientated positive
		$self->{'f2'} = $f  if $f < $self->{'f2'};
		$self->{'t2'} = $t  if $t > $self->{'t2'};
	    } else {
		#label pair orientated negative
		$self->{'f2'} = $f  if $f > $self->{'f2'};
		$self->{'t2'} = $t  if $t < $self->{'t2'};
	    }
	}

	#populate sequence array
	for ($i=0; $i < $len; $i++) {

	    $c = substr($string, $i, 1);

	    #begin/end
	    if ($c eq $Mark_Pad) {
		if ($state == 0) {
		    $self->{'pfxlen'}++;
		} elsif ($state > 0) {
		    $self->{'sfxlen'}++;
		    $state = 2;
		}
		next;
	    }

	    #middle
	    $state = 1    if $state < 1;

	    #skip gaps
	    next    if $c eq $Mark_Gap;

	    $p = $hi - $i;    #REVERSE

	    #warn "insert(-): $p/($self->{'lo'},$self->{'hi'}) = $i/$len\t[$c]\n";

	    #store other text, including Mark_Spc or Mark_Fs[12] symbols
	    if (exists $self->{'seq'}->{$p}) {
		next if $self->{'seq'}->{$p} eq $c;
		warn "${self}::insert(-): assembly clash: [$p] = '$$self->{'seq'}->{$p}' -> '$c'\n"
		    if $Bio::MView::Sequence::WARNCLASH;
		next unless $Bio::MView::Sequence::OVERWRITE;
		#unrecord old frameshift
		$self->{'fs'}--  if $self->{'seq'}->{$p} eq $Mark_Fs1;
		$self->{'fs'}--  if $self->{'seq'}->{$p} eq $Mark_Fs2;
	    }
	    #record frameshift
	    $self->{'fs'}++  if $c eq $Mark_Fs1;
	    $self->{'fs'}++  if $c eq $Mark_Fs2;

	    #write the character
	    $self->{'seq'}->{$p} = $c;
	}

	if ($Bio::MView::Sequence::HARDFASTA and
	    $len != $hi-$lo+1 and $self->{'fs'} < 1) {
	    warn "${self}::insert(-) unframeshifted length mismatch ($len != $hi-$lo+1)\n";
	}
	#warn "insert(-): $self->{'lo'} $self->{'hi'}\n";
    }

    #adjust prefix/suffix positions given new lengths
    $self->{'reflo'} = $self->{'hi'} - $self->{'pfxlen'};
    $self->{'refhi'} = $self->{'lo'} + $self->{'sfxlen'};

    #$self->print;

    $self;
}

sub _substr {
    my ($self, $stop, $start) = @_; #REVERSE

    return ''  if $start < 1 or $stop < 1;   #negative range args
    return ''  unless $self->{'lo'} > 0;     #empty
    return ''  if $stop  > $self->{'hi'};    #missed (too high)
    return ''  if $start < $self->{'lo'};    #missed (too low)

    $start = $self->{'hi'}  if $start > $self->{'hi'};
    $stop  = $self->{'lo'}  if $stop  < $self->{'lo'};
    #warn "_substr(-,lo,hi): $start, $stop";
    #warn "_substr(-,beg,end): $self->{'reflo'}, $self->{'refhi'}";
    $stop--;

    my $s = '';

    for (my $i = $start; $i > $stop; $i--) {
	if ($i < $self->{'refhi'} or $i > $self->{'reflo'}) {
	    $s .= $Mark_Pad;
	    next;
	}
	#warn "_substr(-): $i [$self->{'seq'}->{$i}]";
	if (exists $self->{'seq'}->{$i}) {
	    $s .= $self->{'seq'}->{$i};
	} elsif ($i > $self->{'refhi'} and $i < $self->{'reflo'}) {
	    #We never write a gap at the edge. Frameshifts upset the
	    #numbering so that [begin,end] may exceed the true range.
	    #warn "_substr(-): gap: [$i] ($self->{lo},$self->{hi}), ($self->{reflo},$self->{refhi})\n";
	    $s .= $Mark_Gap;
	}
    }
    $s;
}

sub raw {
    my ($self, $col, $p) = @_;

    $p = $self->{'hi'} - $col + 1;    #REVERSE $col

    return ''  if $p < $self->{'lo'} or $p > $self->{'hi'};

    return $Mark_Pad
	if $p > $self->{'reflo'} or $p < $self->{'refhi'};
    return $self->{'seq'}->{$p}
        if exists $self->{'seq'}->{$p};
    return $Mark_Gap;
}

sub col {
    my ($self, $col, $p) = @_;

    $p = $self->{'hi'} - $col + 1;    #REVERSE $col

#    warn("col(-): $col [", 
#	 exists $self->{'seq'}->{$col} ? $self->{'seq'}->{$col} : '',
#	 "]=> $p [", 
#	 exists $self->{'seq'}->{$p} ? $self->{'seq'}->{$p} : '',
#	 "]\n");

    return ''  if $p < $self->{'lo'} or $p > $self->{'hi'};

    return $self->{'text_pad'}
	if $p > $self->{'reflo'} or $p < $self->{'refhi'};
    if (exists $self->{'seq'}->{$p}) {
	return $self->{'text_spc'}
	    if $self->{'seq'}->{$p} eq $Mark_Spc;
	return $self->{'text_fs1'}
	    if $self->{'seq'}->{$p} eq $Mark_Fs1;
	return $self->{'text_fs2'}
	    if $self->{'seq'}->{$p} eq $Mark_Fs2;
	return $self->{'seq'}->{$p};
    }
    return $self->{'text_gap'};
}


###########################################################################
1;
