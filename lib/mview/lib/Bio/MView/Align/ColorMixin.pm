# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Align::ColorMixin;

use Bio::MView::Color::ColorMap;

use vars qw($FIND_WARNINGS);

$FIND_WARNINGS = 0;           #find block warning count

my $FIND_SEPARATOR = ':';     #for multiple find patterns
my $FIND_COLORMAP  = 'FIND';  #hardwired find colormap

######################################################################
# public methods
######################################################################
sub color_none {
    my ($self, $kw) = @_;

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $kw->{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

        $c = $self->{'string'}->raw($i);

        #warn "[$i]= $c\n";

        #white space: no color
        next  if $self->{'string'}->is_space($c);

        #gap or frameshift: gapcolour
        if ($self->{'string'}->is_non_char($c)) {
            push @$color, $i, 'color' => $kw->{'gapcolor'};
            next;
        }

        push @$color, $i, 'color' => $kw->{'symcolor'};
    }
}

sub color_by_find_block {
    my ($self, $kw) = @_;

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $kw->{'symcolor'};

    my $block = $self->find_blocks($kw->{'find'}, $FIND_COLORMAP);

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

        $c = $self->{'string'}->raw($i);

        #warn "[$i]= $c\n";

        #white space: no color
        next  if $self->{'string'}->is_space($c);

        #gap or frameshift: gapcolour
        if ($self->{'string'}->is_non_char($c)) {
            push @$color, $i, 'color' => $kw->{'gapcolor'};
            next;
        }

        if (exists $block->{$i}) {
            #use symbol color/wildcard colour
            @tmp = $self->get_color($block->{$i}, $FIND_COLORMAP);
        } else {
            @tmp = ();
        }

        push @$color,
            $self->color_tag($kw->{'css1'}, $kw->{'symcolor'}, $i, @tmp);
    }
}

#subclass overrides
sub color_by_type {
    my ($self, $kw) = @_;

    my $color = $self->{'display'}->{'range'};

    push @$color, 1, $self->length, 'color' => $kw->{'symcolor'};

    my $end = $self->length + 1;

    for (my $i=1; $i<$end; $i++) {

        my $c = $self->{'string'}->raw($i);

        #warn "[$i]= $c\n";

        #white space: no color
        next  if $self->{'string'}->is_space($c);

        #gap: gapcolour
        if ($self->{'string'}->is_non_char($c)) {
            push @$color, $i, 'color' => $kw->{'gapcolor'};
            next;
        }

        #use symbol color/wildcard colour
        my @tmp = $self->get_color($c, $kw->{'aln_colormap'});

        push @$color,
            $self->color_tag($kw->{'css1'}, $kw->{'symcolor'}, $i, @tmp);
    }
}

#subclass overrides
sub color_by_identity {
    my $self = shift;
    $self->color_by_identity_body(1, @_);
}

sub color_by_mismatch {
    my $self = shift;
    $self->color_by_identity_body(0, @_);
}

######################################################################
# private methods
######################################################################
sub get_color {
    my ($self, $c, $map) = @_;
    #warn "get_color: $c, $map";

    return $COLORMAP->get_symbol_color($map, $c)
        if $COLORMAP->has_symbol_color($map, $c);

    return $COLORMAP->get_wildcard_color($map)
        if $COLORMAP->has_wildcard_color($map);

    return $COLORMAP->get_palette_color($map, 'S')  #ignore CSS setting
        if $COLORMAP->has_palette_color($map);

    return ();  #no match
}

sub color_tag {
    my ($self, $css, $symcolor, $i) = (shift, shift, shift, shift);
    if (@_) {
        return ($i, 'class' => $_[1] . $_[2])  if $css;
        return ($i, 'color' => $_[0]);
    }
    return ($i, 'color' => $symcolor);
}

sub find_blocks {
    my ($self, $find, $colormap) = @_;

    my $mapsize = $COLORMAP->get_colormap_length($colormap);

    my @patterns = split($FIND_SEPARATOR, $find);

    if (@patterns > $mapsize and $FIND_WARNINGS) {
        warn "find: @{[scalar @patterns]} pattern blocks but only $mapsize color@{[$mapsize gt 1 ? 's' : '']} in colormap '$colormap' - recycling\n";
        $FIND_WARNINGS--;
    }

    my $matches = $self->{'string'}->findall(\@patterns, $mapsize);
    my $index = {};

    foreach my $block (@$matches) {
        $index->{$block->[1]} = $block->[0];
    }

    return $index;
}

sub color_by_identity_body {
    my ($self, $byidentity, $kw, $othr) = @_;

    return  unless defined $othr;

    die "${self}::color_by_identity: length mismatch\n"
        unless $self->length == $othr->length;

    my ($color, $end) = ($self->{'display'}->{'range'}, $self->length+1);

    push @$color, 1, $self->length, 'color' => $kw->{'symcolor'};

    for (my $i=1; $i<$end; $i++) {

        my $c1 = $self->{'string'}->raw($i);
        my $c2 = $othr->{'string'}->raw($i);

        #warn "[$i]= $c1 <=> $c2\n";

        #white space: no color
        next  if $self->{'string'}->is_space($c1);

        #gap or frameshift: gapcolour
        if ($self->{'string'}->is_non_char($c1)) {
            push @$color, $i, 'color' => $kw->{'gapcolor'};
            next;
        }

        my @tmp = ();

        #compare symbols, case-insensitive
        if ($byidentity) {  #identity coloring mode
            @tmp = $self->get_color($c1, $kw->{'aln_colormap'})
                if uc $c1 eq uc $c2; #same symbol or symcolor
        } else {  #mismatch coloring mode
            @tmp = $self->get_color($c1, $kw->{'aln_colormap'})
                if uc $c1 ne uc $c2; #different symbol or symcolor
        }

        push @$color,
            $self->color_tag($kw->{'css1'}, $kw->{'symcolor'}, $i, @tmp);
    }
}

###########################################################################
1;
