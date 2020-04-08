# Copyright (C) 2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Parse::Block;

use Bio::Util::Object;

use vars qw(@ISA);

@ISA = qw(Bio::Util::Object);

my $DEBUG = 0;

sub new {
    my $type = shift;
    my ($key, $start, $stop) = @_;

    my $self = {};
    bless $self, $type;

    $self->{'key'}    = $key;
    $self->{'start'}  = $start;
    $self->{'stop'}   = $stop;
    $self->{'record'} = undef;

    $self;
}

sub str {
    my $self = shift;
    my $s = $self->{'key'};
    $s .= ", " . $self->{'start'};
    $s .= ", " . $self->{'stop'};
    $s .= ", " . "(bytes=" . ($self->{'stop'}-$self->{'start'}) . ")";
    $s .= ", " . (defined $self->{'record'} ? $self->{'record'} : "undef");
    return "{" . $s . "}";
}

sub parse {
    my ($self, $parent, $text) = @_;

    #warn "> Block::parse: ", $self->str(), "\n"  if $DEBUG;

    return $self->{'record'}  if defined $self->{'record'};

    my $class = ref($parent) . "::" . $self->{'key'};

    #warn "block::parse: new $class($parent, $text, $self->{'start'}, $self->{'stop'})\n"
    #    if $DEBUG;

    my $o = $class->new($parent, $text, $self->{'start'}, $self->{'stop'});

    $self->{'record'} = $o;

    #warn "< Block::parse: ", $self->str(), "\n"  if $DEBUG;

    return $o;
}

sub free {
    my $self = shift;

    #warn "> Block::free $self $self->{'key'}\n"  if $DEBUG;

    return  unless defined $self->{'record'};

    $self->{'record'}->free();
    $self->{'record'} = undef;

    #warn "< Block::free $self $self->{'key'}\n"  if $DEBUG;;
}

###########################################################################
1;
