# Copyright (C) 2015 Nigel P. Brown
# $Id: Conservation.pm,v 1.1 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Align::Conservation;

use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Align::Row;

use strict;
use vars qw(@ISA $Debug);

@ISA = qw(Bio::MView::Align::Sequence);

$Debug = 0;

sub new {
    my $type = shift;
    warn "${type}::new() (@_)\n"    if $Debug;
    if (@_ < 1) {
	die "${type}::new() missing arguments\n";
    }
    my ($from, $to, $string) = @_;

    my $self = { %Bio::MView::Align::Sequence::Template };

    $self->{'id'}   = "clustal";
    $self->{'type'} = 'conservation';
    $self->{'from'} = $from;
    $self->{'to'}   = $to;
   
    #encode the new "sequence"
    $self->{'string'} = new Bio::MView::Sequence;
    $self->{'string'}->set_find_pad(' ');
    $self->{'string'}->set_find_pad(' ');
    $self->{'string'}->set_pad(' ');
    $self->{'string'}->set_gap(' ');
    $self->{'string'}->insert([$string, $from, $to]);

    bless $self, $type;

    $self->reset_display;

    $self;
}


###########################################################################
1;
