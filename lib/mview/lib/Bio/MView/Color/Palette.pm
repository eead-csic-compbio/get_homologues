# Copyright (C) 2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::MView::Color::Palette;

sub new {
    my $type = shift;
    my $self = {};
    bless $self, $type;

    $self->{'keys'} = {};
    $self->{'rgb'}  = [];

    $self;
}

######################################################################
# public methods
######################################################################
sub size { return scalar @{$_[0]->{'rgb'}} }

sub insert {
    my ($self, $key, $val) = @_;
    my $last = @{$self->{'rgb'}};

    push @{$self->{'rgb'}}, $val;

    $self->{'keys'}->{$key}  = $last;  #name -> index
    $self->{'keys'}->{$last} = $key;   #index -> name
}

sub has_color {
    my ($self, $key) = @_;
    return 1  if exists $self->{'keys'}->{$key};
    return 0;
}

sub get_index { return $_[0]->_get_key($_[1]) }

sub get_color { return $_[0]->_get_key($_[1]) }

sub get_rgb {
    my ($self, $key) = @_;
    return undef  unless exists $self->{'keys'}->{$key};
    my $i = $self->{'keys'}->{$key};
    return $self->{'rgb'}->[$i];
}

sub get_rgb_at_index {
    my ($self, $index) = @_;
    return $self->{'rgb'}->[$index];
}

sub color_names {
    my $self = shift;
    my @tmp = ();
    for (my $i=0; $i < $self->size; $i++) {
        push @tmp, $self->{'keys'}->{$i};
    }
    return @tmp;
}

######################################################################
# private methods
######################################################################
sub _get_key {
    my ($self, $key) = @_;
    return undef  unless exists $self->{'keys'}->{$key};
    return $self->{'keys'}->{$key};
}

###########################################################################
1;
