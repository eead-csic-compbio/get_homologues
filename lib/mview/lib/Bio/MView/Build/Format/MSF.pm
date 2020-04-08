# Copyright (C) 1997-2015 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::MView::Build::Row::MSF;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 0,   1,    'weight',        'weight',     '5N',       '' ],
    ]
}

sub ignore_columns { ['desc', 'posn1', 'posn2']; }


###########################################################################
package Bio::MView::Build::Format::MSF;

use Bio::MView::Build::Align;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying Bio::Parse::Format parser
sub parser { 'MSF' }

sub parse {
    my $self = shift;
    my ($rank, $use, $id, $wgt, $seq, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    foreach $id (@{$self->{'entry'}->parse(qw(NAME))->{'order'}}) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $id);

        #warn "KEEP: ($rank,$id)\n";

        $wgt = $self->{'entry'}->parse(qw(NAME))->{'seq'}->{$id}->{'weight'};
        $seq = $self->{'entry'}->parse(qw(ALIGNMENT))->{'seq'}->{$id};

        push @hit, new Bio::MView::Build::Row::MSF(
            $rank,
            $id,
            '',
            { 'weight' => $wgt }
        );
        $hit[$#hit]->add_frag($seq);
    }
    #map { $_->dump } @hit;

    #free objects
    $self->{'entry'}->free_parsers(qw(NAME ALIGNMENT));

    return \@hit;
}


###########################################################################
1;
