# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Align::Ruler;

use Bio::MView::Align::Row;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Row);

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    die "${type}::new: missing arguments\n"  if @_ < 2;
    my ($length, $refobj) = @_;

    my $self = new Bio::MView::Align::Row('ruler', '');
    bless $self, $type;

    $self->{'length'} = $length;  #alignment width

    $self->reset_display($refobj);

    $self;
}

######################################################################
# public methods
######################################################################
#override
sub length { return $_[0]->{'length'} }

######################################################################
# private methods
######################################################################
#override
sub reset_display {
    my ($self, $refobj) = @_;

    my @labels = ();

    @labels = $refobj->display_column_labels  if defined $refobj;
    #warn "labels: [@{[join(',', @labels)]}]\n";

    $self->SUPER::reset_display(
        'labels' => \@labels,  #label strings
    );
}

######################################################################
# debug
######################################################################
#overrides Bio::MView::Align::Row::string
sub string { 'ruler length= ' . $_[0]->length }

###########################################################################
1;
