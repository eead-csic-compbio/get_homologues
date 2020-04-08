# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Align::Row;

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    die "${type}::new: missing arguments\n"  if @_ < 2;
    my ($display_type, $uid) = @_;

    my $self = {};
    bless $self, $type;

    $self->{'type'} = $display_type;
    $self->{'uid'}  = $uid ne '' ? $uid : $self;

    $self->reset_display;

    $self;
}

######################################################################
# public methods
######################################################################
sub is_sequence  { 0 }  #subclass overrides
sub is_consensus { 0 }  #subclass overrides
sub is_special   { 0 }  #subclass overrides

sub set_display {
    my $self = shift;
    while (@_) {
        my ($key, $val) = (shift, shift);
        my $ref = ref $val;

        #shallow copy referenced data
        if ($ref eq 'HASH') {
            $self->{'display'}->{$key} = { %$val };
            next;
        }
        if ($ref eq 'ARRAY') {
            $self->{'display'}->{$key} = [ @$val ];
            next;
        }

        #labels hack
        if ($key =~ /^label(\d+)$/) {
            $self->{'display'}->{'labels'}->[$1] = $val;
            next;
        }

        #default
        $self->{'display'}->{$key} = $val;
    }
}

sub uid { return $_[0]->{'uid'} }

sub get_display { return $_[0]->{'display'} }

sub get_label {
    return ''  unless defined $_[0]->{'display'}->{'labels'}->[$_[1]];
    return $_[0]->{'display'}->{'labels'}->[$_[1]];
}

######################################################################
# private methods
######################################################################
#subclass overrides
sub reset_display {
    my $self = shift;
    $self->{'display'} = {
        'type'    => $self->{'type'},
        'range'   => [],  #override: feature column range
        'labels'  => [],  #subclass inserts: empty labels
    };
    $self->set_display(@_);
}

#subclass overrides
sub length { 0 }

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

use Bio::Util::Object qw(dump_self);

sub dump { warn dump_self(@_) }

# subclass overrides
sub string { 'row' }

###########################################################################
1;
