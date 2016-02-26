# -*- perl -*-
# Copyright (C) 1999-2015 Nigel P. Brown
# $Id: BLOCKS.pm,v 1.7 2015/06/14 17:09:04 npb Exp $

###########################################################################
#
# BLOCKS parsing consists of:
#   HEADER
# followed by repeat records of:
#   BLOCK
# interpersed with garbage including a trailing section.
#
###########################################################################
package NPB::Parse::Format::BLOCKS;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);

#delimit full BLOCKS entry
my $BLOCKS_START          = '^BL\d+:\s+.+';
my $BLOCKS_END            = undef;
my $BLOCKS_Null           = '^\s*$';#'

my $BLOCKS_HEADER         = $BLOCKS_START;
my $BLOCKS_HEADERend      = '^Block';
my $BLOCKS_BLOCK          = $BLOCKS_HEADERend;
my $BLOCKS_BLOCKend       = '^\/\/';


#Consume one entry-worth of input on text stream associated with $file and
#return a new BLOCKS instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {

	#start of entry
	if ($line =~ /$BLOCKS_START/o and $offset < 0) {
            $offset = $parent->{'text'}->startofline;
	    next;
	}

	#end of entry
	#if ($line =~ /$BLOCKS_END/o) {
	#    last;
	#}
    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::BLOCKS(undef, $parent->{'text'}, $offset, $bytes);
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

	#HEADER
	if ($line =~ /$BLOCKS_HEADER/o) {
	    $text->scan_until($BLOCKS_HEADERend, 'HEADER');
	    next;
	}

	#BLOCK lines	       	      
	if ($line =~ /$BLOCKS_BLOCK/o) {
	    $text->scan_until($BLOCKS_BLOCKend, 'BLOCK');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next  if $line =~ /$BLOCKS_Null/o;

	#ignore any other line!
	next;

	#terminal line: ignore
	#next  if $line =~ /$BLOCKS_END/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self->test_records(qw(BLOCK));

    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::BLOCKS::HEADER;

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

    while (defined ($line = $text->next_line)) {
	
	#header line
	if ($line =~ /^(BL\d+):\s+(\S+)?/o) {
	    $self->{'ac'} = $1;
	    $self->{'id'} = (defined $2 ? $2 : '');
	}
	
	#default: ignore
	next;
    }

    $self->warn("missing block accession number/identifier")
	unless defined $self->{'ac'};

    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::BLOCKS::BLOCK;

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

    #initialise these:
    $self->{'ac'}    = '';
    $self->{'de'}    = '';
    $self->{'block'} = [];

    while (defined ($line = $text->next_line)) {

	if ($line =~ /^Block\s+(BL\d+\S*)/) {
	    $self->{'ac'} = $1;
	    next;
	}

	if ($line =~ /^ID\s{3}(.*)/) {
	    $self->{'id'} = $1;
	    next;
	}

	if ($line =~ /^AC\s{3}(\S+);\s+(.*)/) {
	    $self->warn("accession number mismatch ($1 cf. $self->{'ac'})")
		unless $self->{'ac'} and $self->{'ac'} eq $1;
	    next;
	}

	if ($line =~ /^DE\s{3}(.*)/) {
	    #allow multiline, just in case
	    $self->{'de'} .= $1;
	    next;
	}

	if ($line =~ /^BL\s{3}(.*)/) {	
	    $self->{'bl'} = $1;
	    next;
	}

	if ($line =~ /^\s*(\S+)\s*\(\s*(\d+)\)\s(\S+)\s+(\d+)/) {
	    push @{$self->{'block'}}, [ $1, $2, $3, $4 ];
	    next;
	}

	#blank line or empty record: ignore
	next  if $line =~ /$BLOCKS_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self->warn("missing block data") if $self->{'ac'} eq '';

    $self->{'de'} = NPB::Parse::Record::strip_english_newlines($self->{'de'});

    $self;#->examine;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n", 'id', $self->{'id'};
    printf "$x%20s -> %s\n", 'ac', $self->{'ac'};
    printf "$x%20s -> %s\n", 'de', $self->{'de'};
    printf "$x%20s -> %s\n", 'bl', $self->{'bl'};
    for (my $i=0; $i < @{$self->{'block'}}; $i++) {
	print "$x$i => [@{$self->{'block'}->[$i]}]\n";
    }
}


###########################################################################
1;
