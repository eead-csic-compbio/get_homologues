# -*- perl -*-
# Copyright (C) 1999-2015 Nigel P. Brown
# $Id: MIPS.pm,v 1.9 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::MIPS;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full MIPS entry
#my $MIPS_START          = '^\s*(?:\S+)?\s+MIPS:';
#my $MIPS_START          = '^(?:PileUp|\s*(?:\S+)?\s+MIPS:)';
my $MIPS_START          = '^>';
my $MIPS_END            = $MIPS_START;

#MIPS record types
my $MIPS_HEADER         = $MIPS_START;
my $MIPS_HEADERend      = '^L;';
my $MIPS_NAME           = $MIPS_HEADERend;
my $MIPS_NAMEend        = '^C;Alignment';
my $MIPS_ALIGNMENT      = '^\s*\d+';
my $MIPS_ALIGNMENTend   = $MIPS_START;
my $MIPS_Null           = '^\s*$';#'


#Consume one entry-worth of input on text stream associated with $file and
#return a new MIPS instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {

	#start of entry
	if ($line =~ /$MIPS_START/o and $offset < 0) {
	    $offset = $parent->{'text'}->startofline;
	    next;
	}

	#consume rest of stream
	last  if $line =~ /$MIPS_END/o;
    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::MIPS(undef, $parent->{'text'}, $offset, $bytes);
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
	if ($line =~ /$MIPS_HEADER/o) {
	    $text->scan_until($MIPS_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#NAME lines	       	      
	if ($line =~ /$MIPS_NAME/o) {    
	    $text->scan_until($MIPS_NAMEend, 'NAME');
	    next;			       	      
	}				       	      
	
	#ALIGNMENT lines		       	      
	if ($line =~ /$MIPS_ALIGNMENT/o) {       	      
	    $text->scan_until($MIPS_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next    if $line =~ /$MIPS_Null/o;

	#end of NAME section: ignore
	next    if $line =~ /$MIPS_NAMEend/o;
	
	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::MIPS::HEADER;

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

    $self->{'desc'} = '';

    #consume Name lines
    while (defined ($line = $text->next_line)) { 

	#> line
	if ($line =~ /^>[^;]+;(\S+)/o) {
	    $self->test_args($line, $1);
	    $self->{'ac'} = $1;
	    next;
	}
	
	#accumulate other lines
	$self->{'desc'} .= $line;
    }

    $self->warn("missing MIPS data\n")  unless exists $self->{'ac'};

    $self->{'desc'} = NPB::Parse::Record::strip_english_newlines($self->{'desc'});

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",   'ac',   $self->{'ac'};
    printf "$x%20s -> '%s'\n", 'desc', $self->{'desc'};
}


###########################################################################
package NPB::Parse::Format::MIPS::NAME;

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

    $self->{'seq'}   = {};
    $self->{'order'} = [];
    
    #consume Name lines
    while (defined ($line = $text->next_line)) {

	if ($line =~ /^L;(\S+)\s+(.*)/o) {
	    $self->test_args($line, $1,$2);
	    $self->{'seq'}->{$1} = NPB::Parse::Record::strip_english_newlines($2);
	    push @{$self->{'order'}}, $1;
	    next;
	} 
	
	next  if $line =~ /$MIPS_Null/;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    foreach my $i (@{$self->{'order'}}) {
	printf "$x%20s -> %-15s %s=%s\n", 
	'seq',  $i,
	'desc', $self->{'seq'}->{$i};
    }
}


###########################################################################
package NPB::Parse::Format::MIPS::ALIGNMENT;

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

    local $^W=0;
    local $_;
    
    $self->{'seq'} = {};

    while (defined ($line = $text->next_line)) {
    
	no strict;

	#start/end positions
	next  if $line =~ /^\s*\d+[^0-9]*\d+\s*$/o;

	#id/sequence
	if ($line =~ /^\s*(\S+)\s+([^0-9]+)\s+\d+$/o) {
	    $self->test_args($line, $1, $2);
	    $self->{'seq'}->{$1} .= $2;
	    next;
	} 

	#default: ignore all other line types (site and consensus data)
    }

    foreach (keys %{$self->{'seq'}}) {
	$self->{'seq'}->{$_} =~ s/ //g;
    }

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    foreach my $i (sort keys %{$self->{'seq'}}) {
	printf "$x%20s -> %-15s =  %s\n", 'seq', $i, $self->{'seq'}->{$i};
    }
}


###########################################################################
1;
