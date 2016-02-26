# -*- perl -*-
# Copyright (C) 1997-2006 Nigel P. Brown
# $Id: Slurp.pm,v 1.15 2005/12/12 20:42:48 brown Exp $

###########################################################################
#
# Used to give Parse behaviour to any unformatted stream
#
###########################################################################
package NPB::Parse::Format::Slurp;

use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#Consume one entry-worth of input on text stream associated with $file and
#return a new Slurp instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    while (defined ($line = $parent->{'text'}->getline)) {

	#start of entry
	if ($offset < 0) {
	    $offset = $parent->{'text'}->startofline;
	    next;
	}

    }
    return 0   if $offset < 0;

    $bytes = $parent->{'text'}->tell - $offset;

    new NPB::Parse::Format::Slurp(undef, $parent->{'text'}, $offset, $bytes);
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
    $self;#->examine;
}


###########################################################################
1;
