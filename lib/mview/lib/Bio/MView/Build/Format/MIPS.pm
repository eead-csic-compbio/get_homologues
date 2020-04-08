# Copyright (C) 1997-2015 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::MView::Build::Format::MIPS;

use Bio::MView::Build::Align;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying Bio::Parse::Format parser
sub parser { 'MIPS' }

sub parse {
    my $self = shift;
    my ($rank, $use, $id, $des, $seq, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    foreach $id (@{$self->{'entry'}->parse(qw(NAME))->{'order'}}) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $id);

        #warn "KEEP: ($rank,$id)\n";

        $des = $self->{'entry'}->parse(qw(NAME))->{'seq'}->{$id};
        $seq = $self->{'entry'}->parse(qw(ALIGNMENT))->{'seq'}->{$id};

        push @hit, new Bio::MView::Build::Row::MIPS($rank, $id, $des, $seq);
    }
    #map { $_->dump } @hit;

    #free objects
    $self->{'entry'}->free_parsers(qw(NAME ALIGNMENT));

    return \@hit;
}


###########################################################################
package Bio::MView::Build::Row::MIPS;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Simple_Row);

sub ignore_columns { ['posn1', 'posn2']; }


###########################################################################
1;
