# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Build::Simple_Row;

use Bio::MView::Build::Row;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

#implements a row with an optional supplied sequence and no info labels
sub new {
    my $type = shift;
    my ($num, $id, $desc, $seq) = @_;

    my $self = new Bio::MView::Build::Row($num, $id, $desc);
    bless $self, $type;

    $self->add_frag($seq)  if defined $seq;

    $self;
}

###########################################################################
package Bio::MView::Build::Row;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::RowInfoMixin;
use Bio::MView::Sequence::Base;
use Bio::MView::SRS;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::RowInfoMixin);

my $DEF_TEXTWIDTH = 30;  #default width to truncate 'text' field

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    die "${type}::new: missing arguments\n"  if @_ < 3;
    my ($num, $id, $desc) = (shift, shift, shift);
    my $self = {};

    bless $self, $type;

    #strip non-identifier leading rubbish: > or /:
    $id =~ s/^(>|\/:)//;

    $self->{'rid'}  = $id;                      #raw identifier
    $self->{'uid'}  = $self->uniqid($num, $id); #unique identifier

    $id =~ s/^\#//;                             #comment: special row?
    $id = ' '  unless length($id) > 0;          #non-null for Build::map_id

    $self->{'cid'}  = $id;                      #cleaned identifier

    $self->{'num'}  = $num;                     #row number/string
    $self->{'desc'} = $desc;                    #description string
    $self->{'frag'} = [];                       #list of fragments

    $self->{'seq'}  = new Bio::MView::Sequence(
        $PAR->get('sequences')                  #finished sequence
    );

    $self->{'url'}  = Bio::MView::SRS::srsLink($self->{'cid'});  #url

    $self->save_info(@_);                       #other info

    $self;
}

######################################################################
# protected methods
######################################################################
#subclass overrides: routine to sort 'frag' list: default is null
sub sort {$_[0]}

######################################################################
# public methods
######################################################################
sub num   { $_[0]->{'num'} }
sub num0  { $_[0]->{'num'} ne '' ? $_[0]->{'num'} : '0' }

sub rid   { $_[0]->{'rid'} }
sub uid   { $_[0]->{'uid'} }
sub cid   { $_[0]->{'cid'} }

sub url   { $_[0]->{'url'} }
sub sob   { $_[0]->{'seq'} }

sub desc  { $_[0]->{'desc'} }  #row description
sub covr  { $_[0]->{'covr'} }  #percent coverage
sub pcid  { $_[0]->{'pcid'} }  #percent identity

sub posn1 { '' }  #subclass overrides: first sequence range
sub posn2 { '' }  #subclass overrides: second sequence range

#return the sequence string, if any
sub seq {
    return ''  unless $PAR->get('sequences');  #no interpolation
    return ''  unless defined $_[0]->{'seq'};
    return $_[0]->{'seq'}->string
}

#return possibly truncated description
sub text {
    my $w = defined $_[1] ? $_[1] : $DEF_TEXTWIDTH;
    $w = length $_[0]->{'desc'}  if $w > length $_[0]->{'desc'};
    return sprintf("%-${w}s", $_[0]->truncate($_[0]->{'desc'}, $w));
}

sub set_coverage { $_[0]->{'covr'} = $_[1] }
sub set_identity { $_[0]->{'pcid'} = $_[1] }

#add a sequence fragment to the 'frag' list with value and positions given
#by first three args. use default positions if called with one arg. other
#optional arguments are special to any subclass of Row.
sub add_frag {
    my $self = shift;
    my ($frag, $qry_from, $qry_to) = (shift, shift, shift);

    $qry_from = 1             unless defined $qry_from;
    $qry_to   = length $frag  unless defined $qry_to;

    push @{$self->{'frag'}}, [ \$frag, $qry_from, $qry_to, @_ ];

    #warn "@{$self->{'frag'}->[-1]}\n";
}

#subclass overrides: compute the maximal positional range of a row
sub range {
    my $self = shift;
    my ($lo, $hi) = ($self->{'frag'}->[0][1], $self->{'frag'}->[0][2]);
    foreach my $frag (@{$self->{'frag'}}) {
        #warn "range: $frag->[1], $frag->[2]\n";
        $lo = $frag->[1]  if $frag->[1] < $lo;
        $lo = $frag->[2]  if $frag->[2] < $lo;
        $hi = $frag->[1]  if $frag->[1] > $hi;
        $hi = $frag->[2]  if $frag->[2] > $hi;
    }
    #warn "range: ($lo, $hi)\n";
    return ($lo, $hi);
}

#convert nucleotide positions to a relative amino acid scale
sub translate_range {
    my ($self, $fm, $to) = @_;
    return (int(($fm+2)/3), int($to/3))   if $fm < $to;  #orientation +
    return (int($fm/3),  int(($to+2)/3))  if $fm > $to;  #orientation -
    die "translate_range: from == to  $fm, $to";
}

#subclass overrides: assemble a row from sequence fragments
sub assemble {
    my ($self, $lo, $hi, $gap) = @_;
    $self->sort;                                 #fragment order
    $self->{'seq'}->insert(@{$self->{'frag'}});  #assemble fragments
    $self->{'seq'}->set_range($lo, $hi);         #set sequence range
    $self->{'seq'}->set_pad($gap);
    $self->{'seq'}->set_gap($gap);
}

#return 1 if supplied id matches own id, 0 otherwise
sub matches_id {
    my ($self, $id) = @_;

    #id is zero indicating query
    if ($id =~ /^0$/) {
        return 1  if $self->{'num'} eq '' or $self->{'num'} eq '0';
        return 0;
    }

    #id is major row number only
    if ($id =~ /^\d+$/) {
        return 1  if $self->{'num'} eq $id;       #match major
        return 1  if $self->{'num'} =~ /^$id\./;  #match major.minor prefix
        return 0;
    }

    #id is major.minor row number pair
    if ($id =~ /^\d+\.\d+$/) {
        return 1  if $self->{'num'} eq $id;
        return 0;
    }

    #id is string identifier
    return 1  if $id eq $self->{'rid'} or $id eq $self->{'cid'};

    #id is regex inside // pair, case-insensitive
    if ($id =~ /^\/.*\/$/) {
        my $r = $id; $r =~ s/^\///; $r =~ s/\/$//;
        return 1  if $self->{'cid'} =~ /$r/i;
        return 0;
    }

    #id is wildcard
    return 1  if $id =~ /^\*$/ or $id =~ /^all$/i;

    #no match
    return 0;
}

######################################################################
# private methods
######################################################################
sub uniqid { "$_[1]\034/$_[2]" }

#truncate a string
sub truncate {
    my ($self, $s, $n) = (@_, $DEF_TEXTWIDTH);
    my $t = substr($s, 0, $n);
    substr($t, -3, 3) = '...'  if length $s > $n;
    return $t;
}

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

use Bio::Util::Object qw(dump_self);

sub dump { warn dump_self(@_) }

###########################################################################
1;
