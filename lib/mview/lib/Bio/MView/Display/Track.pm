# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Display::Track;

use Bio::Util::System qw(vmstat);

sub new {
    my $type = shift;
    my ($parent, $hash) = @_;
    my $self = {};
    bless $self, $type;

    $self->{'parent'}   = $parent;  #parent sequence object

    $self->{'length'}   = $parent->length;
    $self->{'start'}    = $parent->{'start'};
    $self->{'stop'}     = $parent->{'stop'};
    $self->{'forwards'} = $parent->forwards;

    $self->{'string'}   = undef;  #alternative sequence to parent
    $self->{'labels'}   = undef;  #annotation column labels
    $self->{'url'}      = '';     #URL; to associate with label1

    $self->{'string'}   = $hash->{'sequence'} if exists $hash->{'sequence'};
    $self->{'labels'}   = $hash->{'labels'}   if exists $hash->{'labels'};
    $self->{'url'}      = $hash->{'url'}      if exists $hash->{'url'};

    #warn "Track::labels: [@{[join(',', @{$self->{'labels'}})]}]\n";

    $self->{'cursor'} = 0;  #current position in display stream, 1-based

    #are 'parent' and 'string' the same length?
    if (defined $self->{'string'}) {
        my $len = $self->{'string'}->length;
        if ($len != $self->{'length'}) {
            die "${type}::new: parent/sequence length mismatch ($self->{'length'}, $len)\n";
        }
    }

    #character range maps: hold attributes by column
    $self->{'r_color'}  = [];  #colors
    $self->{'r_class'}  = [];  #solid/thin

    $hash->{'range'}    = []  unless exists $hash->{'range'};

    #warn "range: [@{[join(',', @{$hash->{'range'}})]}]\n";

    #populate the range maps
    foreach my $range (@{$self->parse_range($hash->{'range'})}) {
        my $start = $range->{'pos'}->[0];
        my $stop  = $range->{'pos'}->[$#{$range->{'pos'}}] + 1;
        for (my $i = $start; $i < $stop; $i++) {
            foreach my $key (keys %$range) {
                next  if $key eq 'pos';
                $self->{"r_$key"}->[$i] = $range->{$key};
            }
        }
    }

    #vmstat("Track::new: done");
    $self;
}

######################################################################
# public methods
######################################################################
sub is_ruler { 0 }

#return kth label string
sub label {
    return ''  unless defined $_[0]->{'labels'}->[$_[1]];
    return $_[0]->{'labels'}->[$_[1]];
}

#return kth labelwidth
sub labelwidth {
    return 0  unless defined $_[0]->{'labels'}->[$_[1]];
    return length($_[0]->{'labels'}->[$_[1]]);
}

#forward/reverse oriented rulers/sequences
sub reset { $_[0]->{'cursor'} = 1 }

#iterator: subclass overrides
sub next_segment { die "$_[0]::new: virtual function called\n" }

######################################################################
# private methods
######################################################################
sub parse_range {
    my ($self, $data) = @_;
    my ($state, @list, $attr, $i) = (2);

    for ($i=0; $i < @$data; $i++) {

        #position data
        if ($data->[$i] =~ /^\d+$/) {
            if ($state == 2) {
                #warn "SAVE   [$i] $data->[$i]\n";
                #begin new range: save old attribute hash, if useful
                if (defined $attr and keys %$attr > 1) {
                    if ($attr->{'pos'}->[0] > $attr->{'pos'}->[1]) {
                        warn "bad attribute range '@{$attr->{'pos'}}'\n";
                    } else {
                        push @list, $attr;
                    }
                }
                $attr  = {};   #reset attribute hash
                $state = 0;    #ready for first pos
            }
            if ($state == 0) {
                #warn "0->1   [$i] $data->[$i]\n";
                #first position
                $attr->{'pos'}->[0] = $data->[$i];
                $state = 1;    #ready for second pos, if any
                next;
            }
            if ($state == 1) {
                #warn "1->2   [$i] $data->[$i]\n";
                #second position
                $attr->{'pos'}->[1] = $data->[$i];
                $state = 2;
                next;
            }
            die "${self}::parse_range: never get here!\n";
        }

        #first position wasn't given: bad attribute list
        if ($state == 0) {
            warn "$self: bad attribute list '@$data'\n";
            last;
        }

        #second position wasn't given: make one and continue
        if ($state == 1) {
            $attr->{'pos'}->[1] = $attr->{'pos'}->[0];
            $state = 2;
        }

        #position range attributes: key=>value pair
        if (grep /^$data->[$i]$/, qw(color class)) {
            #warn "ATTR   [$i] $data->[$i] $data->[$i+1]\n";
            $attr->{$data->[$i]} = $data->[$i+1]  if defined $data->[$i+1];
            $i++;    #jump value
            next;
        }

        warn "$self: unknown range attribute '$data->[$i]'\n";
        $i++;    #skip likely value
    }

    #save last iteration
    if ($state == 2 and defined $attr) {
        #warn "SAVE   LAST\n";
        #save old attribute hash, if useful
        if (defined $attr and keys %$attr > 1) {
            if ($attr->{'pos'}->[0] > $attr->{'pos'}->[1]) {
                warn "bad attribute range '@{$attr->{'pos'}}'\n";
            } else {
                push @list, $attr;
            }
        }
    }

    #resort ranges by (1) lowest 'start', and (2) highest 'stop'
    return [ sort { my $c = $a->{'pos'}->[0] <=> $b->{'pos'}->[0];
                    return  $c  if $c != 0;
                    return  $b->{'pos'}->[1] <=> $a->{'pos'}->[1];
             } @list
    ];
}

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

use Bio::Util::Object qw(dump_self);

sub dump { warn dump_self(@_) }

###########################################################################
1;
