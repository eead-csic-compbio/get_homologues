# Copyright (C) 2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Align::Special;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Align::Sequence;
use Bio::MView::Color::ColorMap;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Sequence);

######################################################################
# public methods
######################################################################
#override
sub is_sequence { 0 }

#override
sub is_special { 1 }

#override
sub set_display {
    my $self = shift;
    #apply any caller's settings and then silence rownum
    $self->SUPER::set_display(@_, 'label0' => '');
}

#subclass overrides
sub color_special {
    my $self = shift;

    my $kw = $PAR->as_dict;

    my $map = $self->get_special_colormap_for_id($kw, $self->{'uid'});

    return  unless defined $map;

    $kw->{'aln_colormap'} = $map;

    $self->color_by_type($kw);
}

######################################################################
# private methods
######################################################################
#Match row identifier like /#MAP/ or /#MAP:/ or /#MAP:stuff/
#where MAP is some known colormap, and return the colormap;
#note: the internal id may have leading text before the hash.
sub get_special_colormap_for_id {
    my ($self, $kw, $id) = @_;
    my ($size, $map) = (0, undef);
    foreach my $m ($COLORMAP->colormap_names) {
        if ($id =~ /\#$m(|:|:.*)$/i) {
            if (length($&) > $size) {
                $size = length($&);
                $map = $m;
            }
        }
    }
    return $map;
}

###########################################################################
1;
