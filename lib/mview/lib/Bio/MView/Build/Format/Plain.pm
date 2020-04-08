# Copyright (C) 1997-2017 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::MView::Build::Format::Plain;

use Bio::MView::Build::Align;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying Bio::Parse::Format parser
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

        push @hit, new Bio::MView::Build::Row::Plain($rank, $id, '', $seq);
    }
    #map { $_->dump } @hit;

    #free objects
    $self->{'entry'}->free_parsers(qw(ALIGNMENT));

    return \@hit;
}


###########################################################################
package Bio::MView::Build::Row::Plain;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Simple_Row);

sub ignore_columns { ['desc', 'posn1', 'posn2']; }


###########################################################################
1;
