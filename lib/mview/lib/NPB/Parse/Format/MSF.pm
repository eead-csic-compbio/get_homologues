# -*- perl -*-
# Copyright (C) 1998-2015 Nigel P. Brown
# $Id: MSF.pm,v 1.22 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::MSF;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#delimit full MSF entry
#my $MSF_START          = '^\s*(?:\S+)?\s+MSF:';
#my $MSF_START          = '^(?:PileUp|\s*(?:\S+)?\s+MSF:)';
my $MSF_START          = '^(?:PileUp|.*\s+MSF:)';
my $MSF_END            = '^PileUp';

#MSF record types
my $MSF_HEADER         = $MSF_START;
my $MSF_HEADERend      = '^\s*Name:';
my $MSF_NAME           = $MSF_HEADERend;
my $MSF_NAMEend        = '^\/\/';
my $MSF_ALIGNMENT      = '^\s*\S+\s+\S';
my $MSF_ALIGNMENTend   = $MSF_START;
my $MSF_Null           = '^\s*$';#'


#Consume one entry-worth of input on text stream associated with $file and
#return a new MSF instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {

	#start of entry
	if ($line =~ /$MSF_START/o and $offset < 0) {
	    $offset = $parent->{'text'}->startofline;
	    next;
	}

	#consume rest of stream
	last  if $line =~ /$MSF_END/o;
    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::MSF(undef, $parent->{'text'}, $offset, $bytes);
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
	if ($line =~ /$MSF_HEADER/o) {
	    $text->scan_until($MSF_HEADERend, 'HEADER');
	    next;
	}
	
	#consume data

	#NAME lines	       	      
	if ($line =~ /$MSF_NAME/o) {    
	    $text->scan_until($MSF_NAMEend, 'NAME');
	    next;			       	      
	}				       	      
	
	#ALIGNMENT lines		       	      
	if ($line =~ /$MSF_ALIGNMENT/o) {       	      
	    $text->scan_until($MSF_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next    if $line =~ /$MSF_Null/o;

	#end of NAME section: ignore
	next    if $line =~ /$MSF_NAMEend/o;
	
	#default
	$self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::MSF::HEADER;

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

    #consume Name lines
    while (defined ($line = $text->next_line)) { 

	#MSF line
	if ($line =~ /^
	    \s*
	    ((?:.+)?)
	    MSF\:\s+(\d+)
	    \s+
	    Type\:\s+(\S+)
            \s*
	    ((?:.+)?)            
	    Check\:\s+(\d+)
	    \s+\.\.
	    /xo) {

	    $self->test_args($line, $2,$3,$5);
	    (
	     $self->{'file'},
	     $self->{'msf'},
	     $self->{'type'},
             $self->{'data'},
	     $self->{'check'},
	    ) = (NPB::Parse::Record::strip_trailing_space($1),
                 $2,$3,
                 NPB::Parse::Record::strip_trailing_space($4),
                 $5);
	}
	
	#ignore any other text
    }

    $self->warn("missing MSF data\n")  unless exists $self->{'msf'};

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> '%s'\n", 'file',   $self->{'file'};
    printf "$x%20s -> %s\n",   'msf',    $self->{'msf'};
    printf "$x%20s -> %s\n",   'type',   $self->{'type'};
    printf "$x%20s -> '%s'\n", 'data',   $self->{'data'};
    printf "$x%20s -> %s\n",   'check',  $self->{'check'};
}


###########################################################################
package NPB::Parse::Format::MSF::NAME;

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
	my $id = "";
	
	if ($line =~ /^(\s*Name:\s+(.+))Len:\s+\d/) {
	    $line = substr($line, length($1));
	    $id = $2;
	    if ($id =~ /^(.*)\s+oo\s*$/) { #weird clustal insertion
		$id = $1;
	    }
	    $id = NPB::Parse::Record::strip_trailing_space($id);
	    #warn "[$id] [$line]\n";

	    if ($line =~ /^
	        Len\:\s+(\d+)     #sequence length
	        \s+
	        Check\:\s+(\S+)   #checksum
	        \s+
	        Weight\:\s+(\S+)  #sequence weight
	        /xo) {
		
	        $self->test_args($line, $1,$2,$3);
	        $self->{'seq'}->{$id} = {
		    'length' => $1,
		    'check'  => $2,
		    'weight' => $3,
	    	};
	        push @{$self->{'order'}}, $id;
	        next;
	    }
	}  
	
	next  if $line =~ /$MSF_Null/;

	#default
	$self->warn("unknown field: $line");
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    foreach my $i (@{$self->{'order'}}) {
	printf "$x%20s -> %-15s %s=%5s %s=%5s %s=%5s\n", 
	'seq',    $i,
	'length', $self->{'seq'}->{$i}->{'length'},
	'check',  $self->{'seq'}->{$i}->{'check'},
	'weight', $self->{'seq'}->{$i}->{'weight'};
    }
}


###########################################################################
package NPB::Parse::Format::MSF::ALIGNMENT;

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

    #warn "@{[keys %$parent]}";
    my $maxnamelen = 0;
    foreach my $id (@{$parent->{record_by_type}->{NAME}[0][3]->{order}}) {
	my $len = length($id);
	$maxnamelen = $len  if $len > $maxnamelen;
    }
    #warn $maxnamelen;
    
    while (defined ($line = $text->next_line)) {
    
	no strict;

	#start/end positions
	next  if $line =~ /^\s*\d+\s+\d+$/o;

	#end position
	next  if $line =~ /^\s*\d+\s*$/o;

	#id/sequence
	if ($line =~ /^\s*(.{$maxnamelen})\s+(.*)$/o) {
	    $id = NPB::Parse::Record::strip_leading_space($1);
	    $id = NPB::Parse::Record::strip_trailing_space($id);
	    $self->test_args($line, $id, $2);
	    $self->{'seq'}->{$id} .= $2;
	    next;
	} 

	next  if $line =~ /$MSF_Null/;

	#default
	$self->warn("unknown field: $line");
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
