# Copyright (C) 1998-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::MView::Build::Format::JNETZ;

use Bio::MView::Build::Align;
use Bio::MView::Color::ColorMap;
use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying Bio::Parse::Format parser
sub parser { 'JNETZ' }

sub parse {
    my $self = shift;
    my ($rank, $use, $aln, $id, $seq, $row, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    $aln = $self->{'entry'}->parse(qw(ALIGNMENT));

    while (1) {
        $rank++; last  if $self->topn_done($rank);

        $seq = $aln->get_query;
        $row = new Bio::MView::Build::Row::JNETZ($rank, 'res', '', $seq);
        push @hit, $row;

        $rank++; last  if $self->topn_done($rank);

        $seq = $aln->get_align;
        $row = new Bio::MView::Build::Row::JNETZ($rank, 'align', '', $seq);
        push @hit, $row;

        $rank++; last  if $self->topn_done($rank);

        $seq = $aln->get_conf;
        $row = new Bio::MView::Build::Row::JNETZ($rank, 'conf', '', $seq);
        push @hit, $row;

        $rank++; last  if $self->topn_done($rank);

        $seq = $aln->get_final;
        $row = new Bio::MView::Build::Row::JNETZ($rank, 'final', '', $seq);
        push @hit, $row;

        last;
    }

    #free objects
    $self->{'entry'}->free_parsers(qw(ALIGNMENT));

    return \@hit;
}

#override
sub make_alignment {
    shift; return new Bio::MView::Build::Format::JNETZ::Align(@_);
}

#construct a header string describing this alignment
sub header {''}


###########################################################################
package Bio::MView::Build::Row::JNETZ;

use Bio::MView::Build::Row;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Simple_Row);

sub ignore_columns { ['num', 'desc', 'covr', 'pcid', 'posn1', 'posn2']; }


###########################################################################
package Bio::MView::Build::Format::JNETZ::Align;

use Bio::MView::Option::Parameters;  #for $PAR
use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Alignment);

#override
sub make_sequence {
    my ($self, $row) = @_;

    my $cid = $row->cid;
    my $uid = $row->uid;
    my $sob = $row->sob;

    my $subtype = $cid;

    #warn "make_sequence: $cid => @{[$subtype ne '' ? $subtype : \"''\"]}\n";

    return new Bio::MView::Align::JNETZ::Sequence($uid, $sob)
        if $subtype eq 'res';

    return new Bio::MView::Align::JNETZ::Prediction($uid, $sob)
        if $subtype eq 'align';

    return new Bio::MView::Align::JNETZ::Prediction($uid, $sob)
        if  $subtype eq 'final';

    return new Bio::MView::Align::JNETZ::Confidence($uid, $sob)
        if $subtype eq 'conf';

    die "${self}::make_sequence: unknown type '$subtype'\n";
}

#override: change the header text
sub color_header {
    my $self = shift;
    my $s = '';
    $s .= "Residues colored by: property\n";
    $s .= "Structure colored by: type\n";
    return $s;
}

#override: ignore generic colouring schemes: use our own
sub set_color_scheme {
    my $self = shift;
    return  unless $PAR->get('html');
    $self->color_special;
}

#override: propagate colour scheme to row objects
sub color_special {
    my $self = shift;
    foreach my $row (@{$self->{'index2row'}}) {
        next  unless defined $row;
        $row->color_special;
    }
}


###########################################################################
package Bio::MView::Align::JNETZ::Sequence;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Align::Special;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Special);

#override
sub color_special {
    my $self = shift;
    my $kw = $PAR->as_dict('css1' => 0);
    $self->color_by_type($kw);
}


###########################################################################
package Bio::MView::Align::JNETZ::Prediction;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Align::Special;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Special);

#override
sub color_special {
    my $self = shift;
    my $kw = $PAR->as_dict('aln_colormap' => 'JNET.PRED');
    $self->color_by_type($kw);
}


###########################################################################
package Bio::MView::Align::JNETZ::Confidence;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Align::Special;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Special);

#override
sub color_special {
    my $self = shift;
    my $kw = $PAR->as_dict('aln_colormap' => 'JNET.CONF');
    $self->color_by_type($kw);
}


###########################################################################
#the 0 here says don't override any colormap of the same name, to
#allow earler loaded user definitions priority - crude, but it'll do.
Bio::MView::Color::ColorMap::load_colormaps(\*DATA, 0);

1;

__DATA__

[JNET.PRED]
Hh  =>  bright-red      #helix
Ee  =>  bright-blue     #sheet
Ll  =>  dark-green      #coil

[JNET.CONF]
0   ->  gray0           #bad
1   ->  gray1
2   ->  gray2
3   ->  gray4
4   ->  gray6
5   ->  gray8
6   ->  gray10
7   ->  gray12
8   ->  gray14
9   ->  gray15          #good

###########################################################################
