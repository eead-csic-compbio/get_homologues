# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Plain.pm,v 1.12 2005/12/12 20:42:48 brown Exp $

###########################################################################
package Bio::MView::Build::Format::Plain;

use Bio::MView::Build::Align;
use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying NPB::Parse::Format parser
sub parser { 'Plain' }

sub parse {
    my $self = shift;
    my ($rank, $use, $entry, $id, $seq, @hit) = (0);
    
    return  unless defined $self->{scheduler}->next;
    
    $entry = $self->{'entry'};

    foreach $id (@{$entry->parse(qw(ALIGNMENT))->{'id'}}) {

	$rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $id);

	#warn "KEEP: ($rank,$id)\n";

	$seq = $entry->parse(qw(ALIGNMENT))->{'seq'}->{$id};
	
	push @hit, new Bio::MView::Build::Simple_Row($rank, $id, '', $seq);
    }
    #map { $_->print } @hit;

    #free objects
    $self->{'entry'}->free(qw(ALIGNMENT));

    return \@hit;
}


###########################################################################
1;
