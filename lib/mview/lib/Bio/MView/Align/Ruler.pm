# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Ruler.pm,v 1.13 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Align::Ruler;

use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Align::Row;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Row);

sub new {
    my $type = shift;
    my ($length, $refobj) = @_;
    my $self = {};
    bless $self, $type;

    $self->{'type'}   = 'ruler';
    $self->{'length'} = $length;

    $self->reset_display($refobj);

    $self;
}

sub id     { $_[0] }
sub string { '' }
sub length { $_[0]->{'length'} }

sub reset_display {
    my ($self, $refobj) = @_;
    my ($data, $pcid) = ('', '');
    if (defined $refobj) {
	$data = $refobj->head;
	$pcid = $refobj->pcid;
    }
    $self->{'display'} =
	{
	 'type'     => $self->{type},
	 'label0'   => '',
	 'label1'   => '',
	 'label2'   => '',
	 'label3'   => $data,
	 'label4'   => $pcid,
	 'label5'   => '',
	 'label6'   => '',
	 'range'    => [],
	 'number'   => 1,
	};
    $self;
}


###########################################################################
1;
