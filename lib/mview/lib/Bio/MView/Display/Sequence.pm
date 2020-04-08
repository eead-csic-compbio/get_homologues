# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Display::Sequence;

use Bio::MView::Display::Track;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Display::Track);

######################################################################
# public methods
######################################################################
#iterator: subclass overrides
sub next_segment {
    my ($self, $par, $chunksize) = @_;
    #warn "${self}::next_segment\n";

    return undef  if $self->{'cursor'} > $self->{'length'};  #done
    return []     unless $par->{'sequences'};                #sequences off

    #current real position
    my $pos1 = $self->{'start'} + $self->{'cursor'} - 1;

    my $rest  = $self->{'length'} - $self->{'cursor'} + 1;  #length remaining
    my $chunk = ($chunksize < $rest ? $chunksize : $rest);

    #warn "($self->{'length'}, $self->{'cursor'}, $chunk, $rest, ($self->{'start'}, $self->{'stop'}))\n";

    my $string = [];
    my $pos2 = $pos1 - 1;

    for (my $i = 0; $i < $chunk; $i++) {
        $pos2++;

        my $c = $self->char_at($self->{'cursor'});

        push @$string, $par->{'dev'}->process_char($c, $self);

        $self->{'cursor'}++;
    }

    $string = $par->{'dev'}->process_segment($string);

    return [ $pos1, $pos2, join('', @$string) ];
}

######################################################################
# private methods
######################################################################
sub char_at {
    defined $_[0]->{'string'} ? $_[0]->{'string'}->col($_[1]) :
        $_[0]->{'parent'}->{'string'}->col($_[1]);
}

###########################################################################
1;
