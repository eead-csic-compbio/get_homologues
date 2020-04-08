# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Sequence;

use Bio::MView::Sequence::Forward;
use Bio::MView::Sequence::Reverse;

my $Find_Pad = '[-._~]';  #input terminal gap characters
my $Find_Gap = '[-._~]';  #input internal gap characters
my $Find_Spc = '\s';      #input whitespace character
my $Find_Fs1 = '\/';      #input frameshift mark '/' from fast[xy], tfast[axy]
my $Find_Fs2 = '\\\\';    #input frameshift mark '\' from fast[xy], tfast[axy]

use vars qw($Mark_Pad $Mark_Gap $Mark_Spc $Mark_Fs1 $Mark_Fs2);

#used by subclasses
$Mark_Pad     = "\001";   #encoded terminal gap
$Mark_Gap     = "\002";   #encoded internal gap
$Mark_Spc     = "\003";   #encoded whitespace
$Mark_Fs1     = "\004";   #encoded frameshift
$Mark_Fs2     = "\005";   #encoded frameshift

my $Text_Pad  = '.';      #default output terminal gap
my $Text_Gap  = '-';      #default output internal gap
my $Text_Spc  = ' ';      #default output whitespace
my $Text_Fs1  = '/';      #default output frameshift1
my $Text_Fs2  = '\\';     #default output frameshift2

my $WARNCLASH = 0;  #set to 1 to enable reporting of clashing symbols
my $OVERWRITE = 0;  #set to 1 to enable clash overwrite by newest symbol
my $HARDFASTA = 0;  #relax length rules for UV. FASTA input

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
# +query  23456 (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 2,6, 0,0)
# +hit    -23-- (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 3,4, 1,2)
#
# column  12345
# +query  23456 (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 2,6, 0,0)
# -hit    -32-- (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 3,4, 1,2)
#
# Reverse query:
#
# column  12345
# -query  65432 (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 6,2, 0,0)
# +hit    -23-- (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 5,4, 1,2)
#
# column  12345
# -query  65432 (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 6,2, 0,0)
# -hit    -32-- (lo,hi, reflo,refhi, pfxlen,sfxlen) = (2,6, 5,4, 1,2)
#

###########################################################################
sub new {
    my $type = shift;
    my ($save) = (@_, 1);
    #warn "${type}::new(@_)\n";
    my $self = {};

    bless $self, $type;

    $self->{'save'}   = $save;  #flag: save sequences if true
    $self->{'seq'}    = {};     #sequence hash of char

    $self->{'lo'}     = 0;      #absolute start of string
    $self->{'hi'}     = 0;      #absolute end of string

    $self->{'pfxlen'} = 0;      #count of prefix padding
    $self->{'sfxlen'} = 0;      #count of suffix padding
    $self->{'seqlen'} = 0;      #count of sequence chars

    $self->{'reflo'}  = 0;      #first position of sequence data
    $self->{'refhi'}  = 0;      #last  position of sequence data

    $self->{'f1'}     = undef;  #from label of sequence
    $self->{'t1'}     = undef;  #to   label of sequence

    $self->{'f2'}     = undef;  #2nd from label of sequence (optional)
    $self->{'t2'}     = undef;  #2nd to   label of sequence (optional)

    $self->{'fs'}     = 0;      #has frameshifts

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

    $self;
}

sub string {
    my $s = $_[0]->_substr($_[0]->{'lo'}, $_[0]->{'hi'});
    $s =~ s/$Mark_Pad/$_[0]->{'text_pad'}/g;
    $s =~ s/$Mark_Gap/$_[0]->{'text_gap'}/g;
    $s =~ s/$Mark_Spc/$_[0]->{'text_spc'}/g;
    $s =~ s/$Mark_Fs1/$_[0]->{'text_fs1'}/g;
    $s =~ s/$Mark_Fs2/$_[0]->{'text_fs2'}/g;
    return $s;
}

sub sequence {
    my $s = $_[0]->_substr($_[0]->{'lo'}, $_[0]->{'hi'});
    $s =~ s/$Mark_Pad//og;
    $s =~ s/$Mark_Gap//og;
    $s =~ s/$Mark_Spc//og;
    $s =~ s/$Mark_Fs1//og;
    $s =~ s/$Mark_Fs2//og;
    return $s;
}

#return full sequence length; includes special characters: gaps, padding, frameshifts
sub length { $_[0]->lablen }

#return pure sequence length; excludes special characters
sub seqlen { $_[0]->{'seqlen'} }

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
    if ($lo != 0 and $hi != 0) {
        die "$self: range values in wrong order ($lo, $hi)\n"    if $lo > $hi;
        ($self->{'lo'}, $self->{'hi'}) = ($lo, $hi);
    }
}

sub make_forward {
    no strict qw(subs);
    bless $_[0], Bio::MView::Sequence::Forward;
}

sub make_reverse {
    no strict qw(subs);
    bless $_[0], Bio::MView::Sequence::Reverse;
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

#subclass overrides
sub raw { die "$_[0]::raw: virtual function called\n" }

#subclass overrides
sub char_at { die "$_[0]::char_at: virtual function called\n" }

#subclass overrides
sub col { die "$_[0]::col: virtual function called\n" }

#subclass overrides
sub is_forwards { die "$_[0]::is_forwards: virtual function called\n" }

#subclass overrides: called with unknown orientation
sub insert {
    my $self = shift;
    #get direction from first fragment range longer than 1
    my $forward = 1;  #default to forward sequence
    foreach my $frag (@_) {
        $forward = 1, last  if $frag->[1] < $frag->[2];
        $forward = 0, last  if $frag->[1] > $frag->[2];
    }
    if ($forward) {
        $self->make_forward;
    } else {
        $self->make_reverse;
    }
    $self->_insert(@_);
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
    my $ordA = ord('A');

    #warn "findall: @{[$self->string]}\n";
    #warn "findall: blocks=[@{[join(', ', @$blocks)]}]\n";

    my $end = scalar @$blocks;

    for (my $blocknum=0; $blocknum<$end; $blocknum++) {
        #warn "findall: $blocknum, $blocks->[$blocknum]\n";
        my $blockid = chr($ordA + $blocknum % $mapsize);
        my $pattern = $blocks->[$blocknum];

        for (my ($i,$j) = (0,0); $j<$seqlen; $i++, $j++) {
            while ($gapped[$i] eq $gapchar) {
                $i++;  #sync gapped/ungapped string indices
            }
            my $s = substr($sequence, $j);
            if ($s =~ /^($pattern)/i) {  #case insensitive
                my @match = split(//, $1);
                #warn "\nfindall: $i,$j [@match] $s\n";
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

######################################################################
# protected methods
######################################################################
#subclass overrides
sub _insert { die "$_[0]::_insert: virtual function called\n" }

#subclass overrides
sub _substr { die "$_[0]::_substr: virtual function called\n" }

#subclass overrides
sub _get_sequence_index { die "$_[0]::_get_sequence_index: virtual function called\n" }

sub _encode {
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
}

sub _update_sequence_limits {
    my ($self, $lo, $hi) = @_;
    $self->{'lo'} = $lo  if $lo < $self->{'lo'} or $self->{'lo'} == 0;
    $self->{'hi'} = $hi  if $hi > $self->{'hi'} or $self->{'hi'} == 0;
}

sub _update_sequence_labels {
    my ($self, $f, $t, $fs, $ts) = @_;

    if (defined $f and $f != 0) {

        $self->{$fs} = $f  unless defined $self->{$fs};
        $self->{$ts} = $t  unless defined $self->{$ts};

        if ($f < $t or $self->{$fs} < $self->{$ts}) {
            #label pair orientated positive
            $self->{$fs} = $f  if $f < $self->{$fs};
            $self->{$ts} = $t  if $t > $self->{$ts};
        } else {
            #label pair orientated negative
            $self->{$fs} = $f  if $f > $self->{$fs};
            $self->{$ts} = $t  if $t < $self->{$ts};
        }
    }
}

sub _populate_sequence_array {
    my ($self, $o, $state, $len, $lo, $hi, $string) = @_;

    return  unless $self->{'save'};

    $self->_encode($string);

    for (my $i=0; $i < $len; $i++) {

        my $c = substr($$string, $i, 1);

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
        $state = 1  if $state < 1;

        #skip gaps
        next  if $c eq $Mark_Gap;

        my $p = $self->_get_sequence_index($lo, $hi, $i);

        #warn "insert($o): $p/($self->{'lo'},$self->{'hi'}) = $i/$len\t[$c]\n";

        #store other text, including Mark_Spc or Mark_Fs[12] symbols
        if (exists $self->{'seq'}->{$p}) {

            next  if $self->{'seq'}->{$p} eq $c;

            warn "insert($o): assembly clash: [$p] = '$$self->{'seq'}->{$p}' -> '$c'\n"
                if $WARNCLASH;

            next  unless $OVERWRITE;

            #unrecord old frameshift
            $self->{'fs'}--  if $self->{'seq'}->{$p} eq $Mark_Fs1;
            $self->{'fs'}--  if $self->{'seq'}->{$p} eq $Mark_Fs2;
        }
        #record frameshift
        $self->{'fs'}++  if $c eq $Mark_Fs1;
        $self->{'fs'}++  if $c eq $Mark_Fs2;

        #write the character
        $self->{'seq'}->{$p} = $c;
        $self->{'seqlen'}++  if $c ne $Mark_Pad and $c ne $Mark_Gap and
            $c ne $Mark_Spc and $c ne $Mark_Fs1 and $c ne $Mark_Fs2;
    }

    if ($HARDFASTA and $len != $hi-$lo+1 and $self->{'fs'} < 1) {
        warn "insert($o) unframeshifted length mismatch ($len != $hi-$lo+1)\n";
    }
}

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

use Bio::Util::Object qw(dump_self);

sub dump { warn dump_self(@_) }

###########################################################################
1;
