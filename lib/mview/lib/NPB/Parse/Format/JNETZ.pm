# -*- perl -*-
# Copyright (C) 1999-2015 Nigel P. Brown
# $Id: JNETZ.pm,v 1.9 2015/06/14 17:09:04 npb Exp $

###########################################################################
#
# JNET -z format parsing consists of:
#   ALIGNMENT
#
###########################################################################
package NPB::Parse::Format::JNETZ;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);

#delimit full JNETZ entry
my $JNETZ_START          = '^START PRED';
my $JNETZ_END            = '^END PRED';
my $JNETZ_Null           = '^\s*$';#'

my $JNETZ_ALIGNMENT      = $JNETZ_START;
my $JNETZ_ALIGNMENTend   = $JNETZ_END;


#Consume one entry-worth of input on text stream associated with $file and
#return a new JNETZ instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {
	
	#start of entry
	if ($line =~ /$JNETZ_START/o and $offset < 0) {
            $offset = $parent->{'text'}->startofline;
	    next;
	}

	#end of entry
	if ($line =~ /$JNETZ_END/o) {
	    last;
	}
    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::JNETZ(undef, $parent->{'text'}, $offset, $bytes);
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

	#ALIGNMENT lines	       	      
	if ($line =~ /$JNETZ_ALIGNMENT/o) {
	    $text->scan_until($JNETZ_ALIGNMENTend, 'ALIGNMENT');
	    next;			       	      
	}				       	      
	
	#blank line or empty record: ignore
	next  if $line =~ /$JNETZ_Null/o;

	#terminal line: ignore
	next  if $line =~ /$JNETZ_END/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package NPB::Parse::Format::JNETZ::ALIGNMENT;

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
    $self->{'query'} = [];
    $self->{'final'} = [];
    $self->{'conf'}  = [];
    $self->{'align'} = [];

    while (defined ($line = $text->next_line)) {

	if ($line =~
	    /^\s*(\S)
	    \s+(\S)
	    \s+\|
	    \s+(\S)
	    \s+.
	    \s+(\d)
	    /xo
	   ) {
	    push @{ $self->{'query'}}, $1;
	    push @{ $self->{'final'}}, $2;
	    push @{ $self->{'align'}}, $3;
	    push @{ $self->{'conf'}},  $4;
	    next;
	}

	#header line: ignore
	next  if $line =~ /$JNETZ_ALIGNMENT/;

	#blank line or empty record: ignore
	next  if $line =~ /$JNETZ_Null/o;

	#default
	$self->warn("unknown field: $line");
    }

    $self;#->examine;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",   'query', scalar $self->get_query;
    printf "$x%20s -> %s\n",   'final', scalar $self->get_final;
    printf "$x%20s -> %s\n",   'conf',  scalar $self->get_conf;
    printf "$x%20s -> %s\n",   'align', scalar $self->get_align;
}

sub get_query {
    my $self = shift;
    my @a = @{$self->{'query'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_final {
    my $self = shift;
    my @a = @{$self->{'final'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_conf {
    my $self = shift;
    my @a = @{$self->{'conf'}};
    return \@a    if wantarray;
    return join("", @a);
}

sub get_align {
    my $self = shift;
    my @a = @{$self->{'align'}};
    return \@a    if wantarray;
    return join("", @a);
}


###########################################################################
1;
