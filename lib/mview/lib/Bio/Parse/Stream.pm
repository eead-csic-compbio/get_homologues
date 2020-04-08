# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Parse::Stream;

use Bio::Util::Object;
use Bio::Parse::Record;
use Bio::Parse::ReadFile;

use vars qw(@ISA);

@ISA = qw(Bio::Util::Object);

#assumes a stream doesn't mix formats
sub new {
    my $type = shift;
    my ($file, $format) = @_;

    my $self = {};
    bless $self, $type;

    $self->{'file'}   = $file;
    $self->{'format'} = $format;
    $self->{'text'}   = new Bio::Parse::ReadFile($file);
    $self->{'offset'} = $self->{'text'}->get_offset();

    load_parser_class($format);

    $self;
}

sub get_entry {
    my $self = shift;
    #warn "Stream::get_entry: offset= $self->{'offset'}\n";

    $self->{'text'}->reset($self->{'offset'});  #start parsing here

    #warn "Stream: next offset: $self->{'offset'}\n";
    my $parser = new_parser($self->{'text'}, $self->{'format'});

    return undef  unless $parser;

    $self->{'offset'} = $self->{'text'}->get_offset();

    return $parser;
}

sub close { $_[0]->{'text'}->close }

###########################################################################
# private methods / statics
###########################################################################
sub DESTROY { $_[0]->close }

sub load_parser_class {
    my $format = shift;
    ($format = "Bio::Parse::Format::$format") =~ s/::/\//g;
    require "$format.pm";
}

sub new_parser {
    my ($text, $format) = @_;
    my $parser = "Bio::Parse::Format::${format}::get_entry";
    no strict 'refs';
    return $parser->($text);
}

###########################################################################
# debug methods
###########################################################################

# called in the parser test suite
sub print {
    my $self = shift;
    print "Class $self\n";
    foreach my $key (qw(file format)) {
        printf "%16s => %s\n", $key,
            defined $self->{$key} ? $self->{$key} : 'undef';
    }
}

###########################################################################
1;
