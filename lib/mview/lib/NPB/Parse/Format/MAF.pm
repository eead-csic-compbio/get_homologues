# -*- perl -*-
# Copyright (C) 2013-2015 Nigel P. Brown
# $Id: MAF.pm,v 1.2 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::MAF;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);

# See: https://cgwb.nci.nih.gov/FAQ/FAQformat.html

#delimit full MAF entry
my $MAF_START          = '^##maf\s';
my $MAF_END            = $MAF_START;

#MAF record types
my $MAF_Null           = '^\s*$';#'
my $MAF_BLOCK          = '^a\s';
my $MAF_BLOCKend       = $MAF_Null;
my $MAF_HEADER         = $MAF_START;
my $MAF_HEADERend      = "(?:$MAF_BLOCK|$MAF_HEADER)";


#Consume one entry-worth of input on text stream associated with $file and
#return a new MAF instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {

	#start of entry
 	if ($line =~ /$MAF_START/o and $offset < 0) {
	    $offset = $parent->{'text'}->startofline;
	    next;
	}

	#consume rest of stream
        if ($line =~ /$MAF_END/o) {
	    last;
        }
    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::MAF(undef, $parent->{'text'}, $offset, $bytes);
}

#Parse one entry
sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);

    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    while (defined ($line = $text->next_line)) {

	#HEADER lines
	if ($line =~ /$MAF_HEADER/o) {
	    $text->scan_until($MAF_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#BLOCK lines
	if ($line =~ /$MAF_BLOCK/o) {
	    $text->scan_until($MAF_BLOCKend, 'BLOCK');
	    next;
	}

	#blank line or empty record: ignore
	next    if $line =~ /$MAF_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::MAF::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $line, $record);

    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    $self->{'version'} = '';
    $self->{'scoring'} = '';
    $self->{'program'} = '';
    $self->{'options'} = '';

    #consume header lines
    while (defined ($line = $text->next_line)) {

	#first line: always has version: 
	if ($line =~ /^##maf\s+version=(\S+)/o) {
	    $self->{'version'} = $1;
	}
	if ($line =~ /^##maf.*?scoring=(\S+)/o) {
	    $self->{'scoring'} = $1;
	}
	if ($line =~ /^##maf.*?program=(\S+)/o) {
	    $self->{'program'} = $1;
	    next;
	}

	#subsequent comment lines: program parameters
	if ($line =~ /^(#\s+.*)/) {
	    $self->{'options'} .= $1 . "\n";
	}

	#ignore any other text
    }
    chomp($self->{'options'});

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n", 'version', $self->{'version'};
    printf "$x%20s -> %s\n", 'scoring', $self->{'scoring'};
    printf "$x%20s -> %s\n", 'program', $self->{'program'};
    printf "$x%20s -> %s\n", 'options', $self->{'options'};
}


###########################################################################
package NPB::Parse::Format::MAF::BLOCK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    my $type = shift;
    if (@_ < 2) {
	#at least two args, ($offset, $bytes are optional).
	NPB::Message::die($type, "new() invalid arguments (@_)");
    }
    my ($parent, $text, $offset, $bytes) = (@_, -1, -1);
    my ($self, $id, $line, $record);

    $self = new NPB::Parse::Record($type, $parent, $text, $offset, $bytes);
    $text = new NPB::Parse::Record_Stream($self);

    local $^W=0;
    local $_;

    $self->{'row'} = [];

    $self->{'score'} = '';
    $self->{'pass'}  = '';

    my $off = 0;

    while (defined ($line = $text->next_line)) {

	no strict;

	chomp $line;

	#a score=xxxx.yyyy
	if ($line =~ /^a\s+.*?score=(\S+)/o) {
	    $self->{'score'} = $1;
	}

	#a pass=x
	if ($line =~ /^a\s+.*?pass=(\S+)/o) {
	    $self->{'pass'} = $1;
	}

	#s src start size strand srcsize sequence
	if ($line =~ /^s\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)/o) {
	    $self->test_args($line, $1, $2, $3, $4, $5, $6);
	    push @{$self->{'row'}}, {
		'id'      => $1,
		'start'   => $2,
		'size'    => $3,
		'strand'  => $4,
		'srcsize' => $5,
		'seq'     => $6,
	    };
	    $off = length($line) - length($6);
	    next;
	}

	next    if $line =~ /^[aie]\s/o;
	next    if $line =~ /$MAF_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    #line length check
    if (defined $self->{'row'}->[0]) {
	$off = length $self->{'row'}->[0]->{'seq'};
	foreach my $row (@{$self->{'row'}}) {
	    $line = $row->{'seq'};
	    my $len = length $line;
	    #warn "$off, $len, $row->{'id'}\n";
	    if ($len != $off) {
		$self->die("length mismatch for '$row->{'id'}' (expect $off, saw $len):\noffending sequence: [$line]\n");
	    }
	}
    }

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n", 'score', $self->{'score'};
    printf "$x%20s -> %s\n", 'pass',  $self->{'pass'};
    foreach my $row (@{$self->{'row'}}) {
	printf "$x%20s -> %s\n", 'id',      $row->{'id'};
	printf "$x%20s -> %s\n", 'start',   $row->{'start'};
	printf "$x%20s -> %s\n", 'size',    $row->{'size'};
	printf "$x%20s -> %s\n", 'strand',  $row->{'strand'};
	printf "$x%20s -> %s\n", 'srcsize', $row->{'srcsize'};
	printf "$x%20s -> %s\n", 'seq',     $row->{'seq'};
    }
}


###########################################################################
1;
