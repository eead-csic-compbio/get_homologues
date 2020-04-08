# Copyright (C) 2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Parse::BlockKeeper;

use Bio::Util::Object;

use vars qw(@ISA);

@ISA = qw(Bio::Util::Object);

my $DEBUG = 0;

sub new {
    my $type = shift;

    my $self = {};
    bless $self, $type;

    $self->{'record_by_posn'} = [];  #list of records
    $self->{'record_by_type'} = {};  #hash of record lists keyed by type

    $self;
}

# return unordered list of keys
sub keys {
    return keys %{$_[0]->{'record_by_type'}};
}

# test for presence of key
sub has_key {
    return exists $_[0]->{'record_by_type'}->{$_[1]};
}

# get a single block by key and 1-based index
sub get_block {
    my ($self, $key, $num) = (@_, 1);
    return undef  unless exists $self->{'record_by_type'}->{$key};
    return undef  unless defined $self->{'record_by_type'}->{$key}->[$num-1];
    return $self->{'record_by_type'}->{$key}->[$num-1];
}


# get an array of blocks by key
sub get_block_list {
    my ($self, $key) = @_;
    return ()  unless exists $self->{'record_by_type'}->{$key};
    return @{$self->{'record_by_type'}->{$key}};
}

# add new block
sub add {
    my ($self, $block) = @_;
    #warn "> BlockKeeper::add: ", $block->str(), "\n";

    push @{$self->{'record_by_posn'}}, $block;

    my $key = $block->{'key'};

    $self->{'record_by_type'}->{$key} = []
        unless exists $self->{'record_by_type'}->{$key};

    push @{$self->{'record_by_type'}->{$key}}, $block;
}

# remove and return last saved block; free its parsed record if any
sub remove_last {
    my $self = shift;

    return undef  unless @{$self->{'record_by_posn'}};

    my $block = pop @{$self->{'record_by_posn'}};

    pop @{$self->{'record_by_type'}->{$block->{'key'}}};

    $block->free();

    return $block;
}

# free all blocks; set own fields to empty
sub free {
    my $self = shift;
    #warn "> BlockKeeper::free $self\n"  if $DEBUG;

    foreach my $block (@{$self->{'record_by_posn'}}) {
        $block->free();
    }

    $self->{'record_by_posn'} = [];
    $self->{'record_by_type'} = {};

    #warn "< BlockKeeper::free $self\n"  if $DEBUG;
}

# free given blocks leaving data structures intact
sub free_parsers {
    my $self = shift;
    #warn "> BlockKeeper::free_parsers $self\n"  if $DEBUG;

    foreach my $key (@_) {
        next  unless exists $self->{'record_by_type'}->{$key};

        #free the items by type then delete them
        foreach my $block (@{$self->{'record_by_type'}->{$key}}) {
            $block->free();
        }
    }

    #warn "< BlockKeeper::free_parsers $self\n"  if $DEBUG;
}

# return index of give block in its particular slot or -1 if not found
sub lookup_block_number {
    my ($self, $key, $start) = @_;

    my $blocks = $self->{'record_by_type'}->{$key};

    for (my $i=0; $i<@$blocks; $i++) {
        return $i  if $blocks->[$i]->{'start'} == $start;
    }

    return -1;
}

# return number of blocks for given key
sub count {
    my ($self, $key) = @_;

    return 0  unless exists $self->{'record_by_type'}->{$key};

    return scalar @{$self->{'record_by_type'}->{$key}};
}

sub dump_records_by_posn {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x ($indent+2);
    my ($s, %count) = ('');
    foreach my $block (@{$self->{'record_by_posn'}}) {
        my $key = $block->{'key'};
        my $blocks = $self->{'record_by_type'}->{$key};
        my $label = $key;
        if (@$blocks > 1) {
            $label .= '/' . ++$count{$key};
        }
        $s .= sprintf "$x%20s -> [%s]\n", $label,
            join(", ",
                 $block->{'key'},
                 $block->{'start'},
                 $block->{'stop'} - $block->{'start'},
            )
            .
            (defined $block->{'record'} ? ", " . $block->{'record'} : '')
    }
    return $s;
}

###########################################################################
1;
