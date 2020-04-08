# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Display::Ruler;

use Bio::MView::Display::Track;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Display::Track);

sub new {
    my $type = shift;
    my $self = new Bio::MView::Display::Track(@_);
    bless $self, $type;

    if ($self->{'forwards'}) {
        $self->{'charlo'}   = '[';
        $self->{'charhi'}   = ']';
    } else {
        $self->{'charlo'}   = ']';
        $self->{'charhi'}   = '[';
    }

    $self;
}

######################################################################
# public methods
######################################################################
sub is_ruler { 1 }

#overrides
sub next_segment {
    my ($self, $par, $chunksize) = @_;
    #warn "${self}::next_segment\n";

    return undef  if $self->{'cursor'} > $self->{'length'};

    #current real position: by orientation
    my $pos1;
    if ($self->{'forwards'}) {
        $pos1 = $self->{'start'} + $self->{'cursor'} - 1;
    } else {
        $pos1 = $self->{'stop'} - $self->{'cursor'} + 1;
    }

    my $rest  = $self->{'length'} - $self->{'cursor'} + 1;  #length remaining
    my $chunk = ($chunksize < $rest ? $chunksize : $rest);

    #warn "($self->{'length'}, $self->{'cursor'}, $chunk, $rest, ($self->{'start'}, $self->{'stop'}))\n";

    my $string = [];
    my $pos2 = $pos1 + ($self->{'forwards'} ? -1 : +1);

    for (my $i = 0; $i < $chunk; $i++) {

        $self->{'forwards'} ? $pos2++ : $pos2--;

        #ruler ends
        if ($pos2 == $self->{'start'}) {
            push @$string, $self->{'charlo'};
            next;
        }
        if ($pos2 == $self->{'stop'}) {
            push @$string, $self->{'charhi'};
            next;
        }

        #ruler markings
        if ($pos2 % 100 == 0) {
            push @$string, substr($pos2,-3,1);  #third digit from right
            next;
        }
        if ($pos2 % 50 == 0) {
            push @$string, ':';
            next;
        }
        if ($pos2 % 10 == 0) {
            push @$string, '.';
            next;
        }
#       if ($pos2 % 5 == 0) {
#           push @$string, '.';
#           next;
#       }

        push @$string, ' ';
    }
    $self->{'cursor'} += $chunk;

    $string = $par->{'dev'}->process_segment($string);

    return [ $pos1, $pos2, join('', @$string) ];
}

###########################################################################
1;
