# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
# Random access a string in memory assuming any DOS CRLF already converted

package Bio::Parse::ReadText;

use Bio::Util::Object;

use vars qw(@ISA $DEBUG);

@ISA = qw(Bio::Util::Object);

$DEBUG = 0;

sub new {
    my $type = shift;
    #warn "${type}::new:\n"  if $DEBUG;
    Bio::Util::Object::die($type, "new() invalid arguments:", @_)  if @_ < 1;
    my ($text, $base) = (@_, 0);
    my $self = {};
    bless $self, $type;

    $self->{'text'}       = $text;
    $self->{'base'}       = $base;
    $self->{'lastoffset'} = undef;
    $self->{'thisoffset'} = undef;
    $self->{'start'}      = undef;
    $self->{'stop'}       = undef;
    $self->{'extent'}     = length ${$self->{'text'}};

    $self->open();

    $self;
}

sub open {
    my $self = shift;
    #warn "String::open:\n"  if $DEBUG;
    $self->reset($self->{'base'});
}

sub close {
    my $self = shift;
    #warn "String::close:\n"  if $DEBUG;
    $self->{'lastoffset'} = undef;
    $self->{'thisoffset'} = undef;
    $self->{'start'} = undef;
    $self->{'stop'} = undef;
}

sub reset {
    my ($self, $offset) = @_;
    #warn "String::reset($offset)\n"  if $DEBUG;
    $self->{'lastoffset'} = $offset;
    $self->{'thisoffset'} = $offset;
    $self->{'start'} = undef;
    $self->{'stop'} = undef;
}

sub get_offset { $_[0]->{'thisoffset'} }

sub start_count {
    $_[0]->{'start'} = $_[0]->{'lastoffset'};
    $_[0]->{'stop'} = undef;
}

sub stop_count_at_start {  #stop counting at start of line
    $_[0]->{'stop'} = $_[0]->{'lastoffset'};
}

sub stop_count_at_end {    #stop counting at end of line
    $_[0]->{'stop'} = $_[0]->{'thisoffset'};
}

sub get_start {
    return $_[0]->{'start'}  if defined $_[0]->{'start'};
    return $_[0]->{'base'};
}

sub get_stop {
    return $_[0]->{'stop'}  if defined $_[0]->{'stop'};
    return $_[0]->{'thisoffset'};
}

sub getline {
    my $self = shift;
    #warn "String::getline(@_)\n"  if $DEBUG;
    my ($line, $offset) = (@_, $self->{'thisoffset'});

    # validate offset
    if ($offset < 0 or $offset > $self->{'extent'}) {
        die("getline: offset out of range ", $offset);
    } elsif ($offset == $self->{'extent'}) {
        return 0;  #EOF
    }

    my $i = index(${$self->{'text'}}, "\n", $offset);

    if ($i < 0) {
        #read remaining text, if any, until end of string
        $$line = CORE::substr(${$self->{'text'}}, $offset);
    } else {
        #read line upto and including eol; cannot be undef
        $$line = CORE::substr(${$self->{'text'}}, $offset, $i - $offset + 1);
    }

    my $bytes = 0 + length $$line;

    $self->{'lastoffset'} = $offset;
    $self->{'thisoffset'} = $offset + $bytes;

    #warn "String::getline: [$$line]\n"  if $DEBUG;

    return $bytes;
}

sub substr {
    my $self = shift;
    #warn "String::substr(@_)\n"  if $DEBUG;
    my ($offset, $bytes) =
        (@_, $self->{'base'}, $self->{'extent'} - $self->{'base'});

    # validate offset
    if ($offset < 0 or $offset > $self->{'extent'}) {
        die("substr: offset out of range ", $offset);
    } elsif ($offset == $self->{'extent'}) {
        return undef;  #EOF
    }

    # validate bytes; truncate if too long
    if ($bytes < 0) {
        die("substr: bytes out of range ", $bytes);
    } elsif (($offset + $bytes) > $self->{'extent'}) {
        $bytes = $self->{'extent'} - $offset;
    }

    my $buff = CORE::substr(${$self->{'text'}}, $offset, $bytes);

    $bytes = 0 + length $buff;

    $self->{'lastoffset'} = $offset;
    $self->{'thisoffset'} = $offset + $bytes;

    #warn "String::substr: [$buff]\n"  if $DEBUG;

    return $buff;
}

###########################################################################
1;
