# Copyright (C) 1998-2015 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::MView::Build::Format::Pearson;

use Bio::MView::Build::Align;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying Bio::Parse::Format parser
sub parser { 'Pearson' }

sub parse {
    my $self = shift;
    my ($rank, $use, $rec, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    foreach $rec ($self->{'entry'}->parse(qw(SEQ))) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $rec->{'id'});

        #warn "KEEP: ($rank,$rec->{'id'})\n";

        push @hit, new Bio::MView::Build::Row::Pearson($rank,
                                                       $rec->{'id'},
                                                       $rec->{'desc'},
                                                       $rec->{'seq'},
            );
    }
    #map { $_->dump } @hit;

    #free objects
    $self->{'entry'}->free_parsers(qw(SEQ));

    return \@hit;
}


###########################################################################
package Bio::MView::Build::Row::Pearson;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Simple_Row);

sub ignore_columns { ['posn1', 'posn2']; }


###########################################################################
1;
