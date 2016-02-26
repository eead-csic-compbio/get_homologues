# -*- perl -*-
# Copyright (C) 1998-2015 Nigel P. Brown
# $Id: CLUSTAL.pm,v 1.21 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::CLUSTAL;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full CLUSTAL entry
my $CLUSTAL_START          = '^\s*CLUSTAL';
my $CLUSTAL_END            = $CLUSTAL_START;

#CLUSTAL record types
my $CLUSTAL_ALIGNMENT      = '^\s*\S+\s+\S+(?:\s+\d+)?\s*$';
my $CLUSTAL_ALIGNMENTend   = $CLUSTAL_START;
my $CLUSTAL_HEADER         = $CLUSTAL_START;
my $CLUSTAL_HEADERend      = "(?:$CLUSTAL_ALIGNMENT|$CLUSTAL_HEADER)";
my $CLUSTAL_Null           = '^\s*$';#'


#Consume one entry-worth of input on text stream associated with $file and
#return a new CLUSTAL instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {

	#start of entry
 	if ($line =~ /$CLUSTAL_START/o and $offset < 0) {
	    $offset = $parent->{'text'}->startofline;
	    next;
	}

	#consume rest of stream
        if ($line =~ /$CLUSTAL_END/o) {
	    last;
        }
    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::CLUSTAL(undef, $parent->{'text'}, $offset, $bytes);
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
	if ($line =~ /$CLUSTAL_HEADER/o) {
	    $text->scan_until($CLUSTAL_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#ALIGNMENT lines		       	      
	if ($line =~ /$CLUSTAL_ALIGNMENT/o) {
	    $text->scan_until($CLUSTAL_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next    if $line =~ /$CLUSTAL_Null/o;

	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::CLUSTAL::HEADER;

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
    $self->{'major'}   = '';
    $self->{'minor'}   = '';

    #consume Name lines
    while (defined ($line = $text->next_line)) { 

	#first part of CLUSTAL line
	if ($line =~ /^
	    \s*
	    CLUSTAL
	    \s+
	    (([^\(\s]+)    #major version, eg., W
	    \s*
	    \((\S+)\))     #minor version, eg., 1.70
	    /xo) {

	    $self->test_args($line, $1, $2, $3);
	    (
	     $self->{'version'},
	     $self->{'major'},
	     $self->{'minor'},
	    ) = ($1, $2, $3);
	    
	}
	
	#first part of T-COFFEE line
	if ($line =~ /^
	    \s*
	    CLUSTAL\sFORMAT\sfor\s((T-COFFEE)
	    \s+
	    Version_(\S+)),     #version number, eg., 1.37
	    /xo) {

	    $self->test_args($line, $1, $2, $3);
	    (
	     $self->{'version'},
	     $self->{'major'},
	     $self->{'minor'},
	    ) = ($1, $2, $3);
	    
	}
	
	#ignore any other text
    }

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",   'version', $self->{'version'};
    printf "$x%20s -> %s\n",   'major',   $self->{'major'};
    printf "$x%20s -> %s\n",   'minor',   $self->{'minor'};
}


###########################################################################
package NPB::Parse::Format::CLUSTAL::ALIGNMENT;

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
    
    $self->{'id'}    = [];
    $self->{'seq'}   = {};
    $self->{'match'} = '';

    my $off = 0;

    while (defined ($line = $text->next_line)) {
    
	no strict;

	chomp $line;

	#match symbols, but only if expected
	if ($off and $line !~ /[^*:. ]/) {
	    $line = substr($line, $off);
	    $self->{'match'} .= $line;
	    $off = 0;
	    next;
	}

	#id/sequence (/optional number)
	if ($line =~ /^\s*(\S+)\s+(\S+)(?:\s+\d+)?$/o) {
	    $self->test_args($line, $1, $2);
	    push @{$self->{'id'}}, $1    unless exists $self->{'seq'}->{$1};
	    $self->{'seq'}->{$1} .= $2;
	    $off = length($line) - length($2);
	    next;
	} 

	next    if $line =~ /$CLUSTAL_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    #line length check (ignore 'match' as this may be missing)
    if (defined $self->{'id'}->[0]) {
	$off = length $self->{'seq'}->{$self->{'id'}->[0]};
	foreach $id (keys %{$self->{'seq'}}) {
	    $line = $self->{'seq'}->{$id};
	    my $len = length $line;
	    #warn "$off, $len, $id\n";
	    if ($len != $off) {
		$self->die("length mismatch for '$id' (expect $off, saw $len):\noffending sequence: [$line]\n");
	    }
	}
    }

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    foreach my $i (@{$self->{'id'}}) {
	printf "$x%20s -> %-15s %s\n", 'seq', $i, $self->{'seq'}->{$i};
    }
    printf "$x%20s -> %-15s %s\n", 'match', '', $self->{'match'};
}


###########################################################################
1;
