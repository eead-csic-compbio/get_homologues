# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
# Random access a file converting any DOS CRLF line endings on the fly

package Bio::Parse::ReadFile;

use FileHandle;
use POSIX;  # for SEEK_SET

use Bio::Util::Object;

use vars qw(@ISA $DEBUG);

@ISA = qw(Bio::Util::Object);

$DEBUG = 0;

use vars qw($OPEN $CLOSED $ERROR);

$OPEN = 1; $CLOSED = 2; $ERROR = 4;

sub new {
    my $type = shift;
    #warn "${type}::new\n"  if $DEBUG;
    Bio::Util::Object::die($type, "new() invalid arguments:", @_)  if @_ < 1;
    my ($file, $base) = (@_, 0);
    my $self = {};
    bless $self, $type;

    $self->{'file'}       = $file;
    $self->{'base'}       = $base;
    $self->{'lastoffset'} = undef;
    $self->{'thisoffset'} = undef;
    $self->{'start'}      = undef;
    $self->{'stop'}       = undef;
    $self->{'extent'}     = undef;

    $self->{'fh'}    = $self->{'fh'} = new FileHandle();
    $self->{'state'} = $CLOSED;
    $self->{'dos'}   = 0;  #true if file is non-native DOS

    $self->open();

    $self;
}

sub open {
    my $self = shift;
    #warn "File::open:\n"  if $DEBUG;
    if ($self->{'state'} != $CLOSED) {
        $self->close();
    }
    $self->{'fh'}->open($self->{'file'}) or $self->die("open: can't open '$self->{'file'}'");
    $self->{'extent'} = (stat $self->{'fh'})[7];
    $self->{'state'} = $OPEN;
    $self->{'dos'} = $self->_file_has_crlf();

    $self->reset($self->{'base'});
}

sub close {
    my $self = shift;
    #warn "File::close:\n"  if $DEBUG;
    if ($self->{'state'} == $CLOSED) {
        return;
    }
    $self->{'fh'}->close();
    $self->{'state'} = $CLOSED;
    $self->{'lastoffset'} = undef;
    $self->{'thisoffset'} = undef;
    $self->{'start'} = undef;
    $self->{'stop'} = undef;
}

sub reset {
    my ($self, $offset) = @_;
    #warn "File::reset($offset)\n"  if $DEBUG;
    if ($self->{'state'} != $OPEN) {
        $self->die("reset: can't reset '$self->{'file'}'");
    }
    $self->{'fh'}->seek($offset, SEEK_SET);  #absolute
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

    if ($self->{'state'} != $OPEN) {
        $self->die("getline: can't read '$self->{'file'}'");
    }

    my ($line, $offset) = (@_, $self->{'thisoffset'});

    #warn sprintf("File::getline: > lastoffset=%d thisoffset=%d offset=%d\n",
    #             $self->{'lastoffset'}, $self->{'thisoffset'}, $offset)  if $DEBUG;

    # validate offset
    if ($offset < 0 or $offset > $self->{'extent'}) {
        die("File::getline: offset ", $offset, " beyond extent 0..", $self->{'extent'});
    } elsif ($offset == $self->{'extent'}) {
        return 0;  #EOF
    }

    my $fh = $self->{'fh'};

    #seek?
    if ((my $delta = $offset - $self->{'thisoffset'}) != 0) {
        #warn "File::getline:  seek($delta)\n"  if $DEBUG;
        #seek relative character counts are wrong on DOS
        return 0  unless seek($fh, $offset, SEEK_SET);  #absolute
    }

    #warn "File::getline:   reading from core::readline()\n";
    $$line = readline($fh);

    return 0  unless defined $$line;

    my $bytes = length $$line;

    $self->{'lastoffset'} = $offset;
    $self->{'thisoffset'} = $offset + $bytes;

    #convert CRLF to LF if file from DOS
    CORE::substr($$line, $bytes-2, 1, '')  if $self->{'dos'};

    #warn "File::getline:  [$$line]\n"  if $DEBUG;

    #warn sprintf("File::getline: < lastoffset=%d thisoffset=%d bytes=%d buff=[%s]\n", $self->{'lastoffset'},
    #             $self->{'thisoffset'}, $bytes, $$line)  if $DEBUG;

    return $bytes;
}

sub substr {
    my $self = shift;
    #warn "File::substr(@_)\n"  if $DEBUG;

    if ($self->{'state'} != $OPEN) {
        $self->die("substr: can't read '$self->{'file'}'");
    }

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

    my $fh = $self->{'fh'};

    if ((my $delta = $offset - $self->{'thisoffset'}) != 0) {
        #warn "File::substr:  seek($delta)\n"  if $DEBUG;
        #seek relative character counts are wrong on DOS
        return undef  unless seek($fh, $offset, SEEK_SET);  #absolute
    }

    my $buff = ''; $bytes = read($fh, $buff, $bytes);

    $self->{'lastoffset'} = $offset;
    $self->{'thisoffset'} = $offset + $bytes;

    #strip multiple CR if file from DOS
    $buff =~ s/\015\012/\012/go  if $self->{'dos'};

    #warn "File::substr: [$buff]\n"  if $DEBUG;

    return $buff;
}

###########################################################################
# private methods
###########################################################################
#detect CRLF in non-native DOS file on UNIX
sub _file_has_crlf {
    my $self = shift;
    my $line = readline($self->{'fh'});
    my $test = index($line, "\r\n") > -1;
    #warn "File::_file_has_crlf: [$line] --> @{[$test > 0 ? 'yes' : 'no']}\n";
    return $test;
}

###########################################################################
1;
