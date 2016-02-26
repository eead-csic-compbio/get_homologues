# -*- perl -*-
# Copyright (C) 1996-2006 Nigel P. Brown
# $Id: Index.pm,v 1.32 2005/12/12 20:42:48 brown Exp $

###########################################################################
package NPB::Parse::Index;

use vars qw(@ISA $VERBOSE);

use Fcntl;
use NDBM_File;
use FileHandle;
use NPB::Parse::Message;
use strict;

@ISA = qw(NPB::Parse::Message);

#public using fully qualified package name
$VERBOSE = 0;

#private
my $EXT = 'rec';

#static method: can we read existing database files?
sub readable {
    if (@_ < 1) {
	NPB::Parse::Message::die('NPB::Parse::Index',
				   "readable() missing filename prefix");
    }
    if (-r "$_[0].dir" and -r "$_[0].pag" and -r "$_[0].$EXT") {
	return 1;
    }
    return 0;
}

#static method: can we write existing database files?
sub writeable {
    if (@_ < 1) {
	NPB::Parse::Message::die('NPB::Parse::Index',
				   "writeable() missing filename prefix");
    }
    if (-w "$_[0].dir" and -w "$_[0].pag" and -w "$_[0].$EXT") {
	return 1;
    }
    return 0;
}

#Associate a tied hash with a list of raw database files through an
#ndbmfile prefixed by $dbm. Three calling modes, $rwa = { 'r', 'w', 'a' }
#to {read,create,append}. First three args required ($class, $dbm, $rwa)
#while subsequent arg ($force) optional. $class specifies the underlying
#Class of database entries to use, $force > 1 permits overwrite of existing
#dbm files.
sub new {
    my $type = shift;
    my $self = {};
    bless $self, $type;

    if (@_ < 3) {
	$self->die("new() invalid argument list (@_)");
    }
    my ($class, $dbm, $rwa, $force) = @_;
    my ($file, $flags, $mode, $entry);

    if ($rwa eq 'w') {
	#check for the existence of any existing file for writing
	if (-e "$dbm.dir" or -e "$dbm.pag" or -e "$dbm.$EXT") {
	    if (! $force) {
		$self->warn("new() Prefix $dbm libraries already exist:");
		system("/bin/ls -C $dbm.dir $dbm.pag $dbm.$EXT 1>&2");
		exit 1;
	    }
	}
	if ($force) {
	    unlink("$dbm.dir", "$dbm.pag", "$dbm.$EXT");
	    if (-e "$dbm.dir" or -e "$dbm.pag" or -e "$dbm.$EXT") {
		$self->die("new() can't unlink $dbm.{pag,dir,rec}");
	    }
	}
	($flags, $mode) = (O_CREAT|O_WRONLY, 0640);
	
    } elsif ($rwa eq 'a') {
	#check for the existence of all input files for appending
	if (-w "$dbm.dir" and -w "$dbm.pag" and -w "$dbm.$EXT") {
	    ($flags, $mode) = (O_WRONLY, 0640);
	} else {
	    ($flags, $mode) = (O_CREAT|O_WRONLY, 0640);
	}
	
    } elsif ($rwa eq 'r') {
	#check for the existence of all input files for reading
	if (! (-r "$dbm.dir" and -r "$dbm.pag" and -r "$dbm.$EXT")) {
	    $self->die("new() Prefix $dbm libraries missing!");
	}
	($flags, $mode) = (O_RDONLY, 0440);
    }
	
    #make a dbmfile tied hash
    #O_RDWR capability fails to allow read on SGI, 'IRIX Release 5.3 IP22'
    if (! tie(%{$self->{'hash'}}, 'NDBM_File', $dbm, $flags, $mode)){ 
	$self->die("new() can't open $dbm.{pag,dir}");
    }

    $self->{'class'}      = 'NPB::Parse::' . $class;
    $self->{'dbm'}        = $dbm;
    $self->{'mode'}       = $rwa;
    $self->{'index_file'} = "$dbm.$EXT";
    $self->{'index_fh'}   = new FileHandle;
    $self->{'data_fh'}    = {};

    if ($rwa eq 'r') {
	$self->{'index_fh'}->open("<  $self->{'index_file'}");
    } elsif ($rwa eq 'w') {
	$self->{'index_fh'}->open(">  $self->{'index_file'}");
    } else {
	#assume append
	$self->{'index_fh'}->open(">> $self->{'index_file'}");
    }
    if (! defined $self->{'index_fh'}) {
	$self->die("new() can't open dbm record file $self->{'index_file'}");
    }
    warn "NPB::Parse::Index opened dbm  file ($self->{'class'}, $self->{'dbm'}, '$self->{'mode'}')\n"
	if $VERBOSE;
    $self;
}

sub DESTROY {
    my $self = shift;

    #close the filehandle onto the indices file
    if (defined $self->{'index_fh'}) {
	$self->{'index_fh'}->close;
	undef $self->{'index_fh'};
	warn "NPB::Parse::Index closed dbm  file ($self->{'class'}, $self->{'dbm'}, '$self->{'mode'}')\n"
	    if $VERBOSE;
    }
    
    #close any filehandles onto the raw data files
    my $fh; foreach $fh (keys %{$self->{'data_fh'}}) {
	if (defined $self->{'data_fh'}->{$fh}) {
	    $self->{'data_fh'}->{$fh}->close;
	    undef $self->{'data_fh'}->{$fh};
	    warn "NPB::Parse::Index closed data file ($fh, 'r')\n"
		if $VERBOSE;
	}
    }
}

sub close {
    $_[0]->DESTROY;
}

#Store 2 ints in dbmfile: offset,bytecount in data file.
sub _write_hash {
    my ($self, $entry) = @_;
    my ($key, $indirect, $val);

    $val = $entry->pack_hash();
    $indirect = pack('L2', $self->{'index_fh'}->tell, length($val));
    foreach $key ($entry->get_indices) {
	$self->{'hash'}->{$key} = $indirect;
    }
    $self->{'index_fh'}->print($val);
}

sub each {
    each %{$_[0]->{'hash'}};
}

sub exists {
    exists $_[0]->{'hash'}->{$_[1]};
}

sub keys {
    keys %{$_[0]->{'hash'}};
}

sub values {
    values %{$_[0]->{'hash'}};
}

sub put_entry {
    my $self = shift;
    if (@_ != 1) {
	$self->die("put_entry() invalid argument list (@_)");
    }
    my ($entry) = @_;
    $self->_write_hash($entry);
}

sub get_entry {
    my $self = shift;
    if (@_ != 1) {
	$self->die("get_entry() invalid argument list (@_)");
    }
    my $key = shift;
    my ($val, $offset, $bytes, $entry) = ('');

    ($offset, $bytes) = unpack('L2', $self->{'hash'}->{$key});

    #now seek
    if ($self->{'index_fh'}->seek($offset, 0) < 1) {
	$self->die("get_entry() can't match $key in $self->{'index_file'}");
    }

    #now read
    if (! defined read($self->{'index_fh'}, $val, $bytes)) {
	$self->die("get_entry() can't read() in $self->{'index_file'}");
    }

    $entry = NPB::Parse::Entry->new($self->{'class'});

#    $entry = "new NPB::Parse::Entry(\$self->{'class'})";
#    $entry = eval $entry;
#    if ($@) {
#	chomp $@;
#	$self->die("get_entry() eval failed with error:  [$@]");
#    }

    $entry->unpack_hash($val);
    $self->_extract_entry($entry);
    $entry;
}

sub _extract_entry {
    my ($self, $entry) = @_;
    my ($file, $fh, $text) = ($entry->get_file, 0, 0);

    #open a new stream if one not already open on this raw data file
    if (exists($self->{'data_fh'}->{$file})) {
	$fh = $self->{'data_fh'}->{$file};
    } else {	
	$fh = new FileHandle;

	if ($fh->open("< $file") != 1) {
	    $self->die("get_entry() can't open datafile $file");
	}
	warn "NPB::Parse::Index opened data file ($file, 'r')\n" if $VERBOSE;
	$self->{'data_fh'}->{$file} = $fh;
    }

    #now seek
    if ($fh->seek($entry->get_offset, 0) < 1) {
	$self->die("get_entry() can't seek() in $file");
    }

    #and read
    if (! defined read($fh, $text, $entry->get_bytes)) {
	$self->die("get_entry() can't read() in $file");
    }

    $entry->set_text(\$text);
}


###########################################################################
1;
