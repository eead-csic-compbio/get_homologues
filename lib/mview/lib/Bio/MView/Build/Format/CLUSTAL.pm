# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: CLUSTAL.pm,v 1.11 2005/12/12 20:42:48 brown Exp $

###########################################################################
package Bio::MView::Build::Format::CLUSTAL;

use Bio::MView::Build::Align;
use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying NPB::Parse::Format parser
sub parser { 'CLUSTAL' }

sub parse {
    my $self = shift;
    my ($rank, $use, $id, $seq, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    foreach $id (@{$self->{'entry'}->parse(qw(ALIGNMENT))->{'id'}}) {

	$rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $id);

	#warn "KEEP: ($rank,$id)\n";

	$seq = $self->{'entry'}->parse(qw(ALIGNMENT))->{'seq'}->{$id};

	push @hit, new Bio::MView::Build::Simple_Row($rank, $id, '', $seq);
    }
    #map { $_->print } @hit;

    #free objects
    $self->{'entry'}->free(qw(ALIGNMENT));

    return \@hit;
}


###########################################################################
1;
