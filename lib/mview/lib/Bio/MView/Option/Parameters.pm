# Copyright (C) 2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Option::Parameters;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw($PAR);

use vars qw($PAR);

$PAR = undef;  #unique global instance

sub new {
    my ($type, $par) = @_;
    if (defined $PAR) {
        die "Bio::MView::Option::Parameters instance already exists\n";
    }
    my $self = {};
    bless $self, $type;

    $self->{'p'} = $par;

    return $PAR = $self;
}

######################################################################
# public methods
######################################################################
#Return internal hash cloned and modified with any extra (key,val) pairs
#supplied by the caller.
sub as_dict {
    my $self = shift;
    my %clone = %{$self->{'p'}};
    while (@_) {
        my ($k, $v) = (shift, shift);
        $clone{$k} = $v;
    }
    return \%clone;
}

#Return value at $key or undef if it doesn't exist.
sub get {
    my ($self, $key) = @_;
    return $self->{'p'}->{$key}  if exists $self->{'p'}->{$key};
    return undef;  #no match
}

#Set value at $key to $val and return previous value or undef.
sub set {
    my ($self, $key, $val) = @_;
    my $old = undef;
    $old = $self->{'p'}->{$key}  if exists $self->{'p'}->{$key};
    $self->{'p'}->{$key} = $val;  #set new
    return $old;  #return old
}

#Copy value at key $src to key $dst. Do nothing if key $src doesn't exist.
#Return original value of key $dst or undef.
sub get_set {
    my ($self, $src, $dst) = @_;
    my $old = undef;
    $old = $self->{'p'}->{$dst}  if exists $self->{'p'}->{$dst};
    if (exists $self->{'p'}->{$src}) {
        $self->{'p'}->{$dst} = $self->{'p'}->{$src};
    } else {
        ;  #src doesn't exist: do nothing
    }
    return $old;
}

######################################################################
# debug
######################################################################
use Bio::Util::Object qw(dump_hash);

#Return pretty-printed parameter listing; operate on supplied keys or
#all keys if none given.
sub dump { my $self = shift; return dump_hash($self->{'p'}, @_) }

###########################################################################
1;
