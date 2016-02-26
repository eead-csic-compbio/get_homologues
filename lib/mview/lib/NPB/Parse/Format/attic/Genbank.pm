# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: Genbank.pm,v 1.49 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::Genbank;

use NPB::Parse::FT;
use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);;


#Exported constants
#skip these many chars at front of raw records
$NPB::Parse::Format::Genbank::skip          = 12;

#indents for nested records
$NPB::Parse::Format::Genbank::indent        = 2;

#special Genbank record fields
$NPB::Parse::Format::Genbank::id_entryname  = '[A-Z]\w{1,8}';    #### CHANGE?????? 
$NPB::Parse::Format::Genbank::id_date       = '[0-9]{2}-[A-Z]{3}-[0-9]{4}';
$NPB::Parse::Format::Genbank::ac_accession  = $NPB::Parse::FT::ac_accession;


#Genbank record types
my $Genbank_LOCUS         = '^LOCUS';       
my $Genbank_ACCESSION     = '^ACCESSION';   
my $Genbank_DEFINITION    = '^DEFINITION';  
my $Genbank_KEYWORDS      = '^KEYWORDS';    
my $Genbank_SOURCE        = '^SOURCE';      
my $Genbank_REFERENCE     = '^REFERENCE';   
my $Genbank_COMMENT       = '^COMMENT';     
my $Genbank_FEATURES      = '^FEATURES';    
my $Genbank_BASECOUNT     = '^BASE\ COUNT';   
my $Genbank_End           = '^\/\/';        # end of entry: two slashes
my $Genbank_Null          = '^\S*\s*$';#'   # blank line or empty record


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new Genbank instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {
	
	#LOCUS line: start of entry
	if ($line =~ /$Genbank_LOCUS/o and $offset < 0) {
	    $offset = $fh->tell - length($line);
	    next;
	}

	#// line: end of entry
	if ($line =~ /$Genbank_End/o) {
	    last;
	}
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::Genbank(undef, $text, $offset, $bytes);
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

    my $ob;

    while (defined ($line = $text->next_line)) {

 	#LOCUS line
 	if ($line =~ /$Genbank_LOCUS/o) {
 	    $record = $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'LOCUS');
 	    $ob = new NPB::Parse::Format::Genbank::LOCUS($self, \$record);
 	    $self->push_indices($ob->{'id'});
 	    next;
 	}
	
 	#ACCESSION line
 	if ($line =~ /$Genbank_ACCESSION/o) {
 	    $record = $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'ACCESSION');
 	    $ob = new NPB::Parse::Format::Genbank::ACCESSION($self, \$record);
	    $self->set_relative_key($ob->{'ac'});
	    ##$self->push_indices(@{$ob->{'accessions'}});
 	    next;
 	}
	
 	#DEFINITION line
 	if ($line =~ /$Genbank_DEFINITION/o) {
 	    $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'DEFINITION');
 	    next;
 	}
	
 	#KEYWORDS line
 	if ($line =~ /$Genbank_KEYWORDS/o) {
 	    $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'KEYWORDS');
 	    next;
 	}
	
 	#SOURCE line
 	if ($line =~ /$Genbank_SOURCE/o) {
 	    $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'SOURCE');
 	    next;
 	}
	
 	#REFERENCE line
 	if ($line =~ /$Genbank_REFERENCE/o) {
 	    $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'REFERENCE');
 	    next;
 	}
	
 	#COMMENT line
 	if ($line =~ /$Genbank_COMMENT/o) {
 	    $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'COMMENT');
 	    next;
 	}
	
 	#FEATURES line
 	if ($line =~ /$Genbank_FEATURES/o) {
 	    $text->scan_nest($NPB::Parse::Format::Genbank::indent, 'FEATURES');
 	    next;
 	}
	
 	#BASECOUNT line
 	if ($line =~ /$Genbank_BASECOUNT/o) {
 	    $text->scan_until($Genbank_End, 'BASECOUNT');
	    last;
 	}
	
	#blank line or empty record: ignore
	if ($line =~ /$Genbank_Null/o) {
	    next;
	}

	#default
	$self->warn("unknown field: ", substr($line, 0, $NPB::Parse::Format::Genbank::skip));
    }
    $self;
}


###########################################################################
package NPB::Parse::Format::Genbank::LOCUS;

use vars qw(@ISA);
use NPB::Parse::Regexps;

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
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::Genbank::skip);

    $line = $text->next_line;
    
    if ($line =~ /^\s*
        ($NPB::Parse::Format::Genbank::id_entryname)                        #id
        \s*
        ($RX_Uint)                                               #length
        \s*
        ([^.]+)                                                  #units
        \s*
        (ds-DNA)                                                 #molecule
        \s*
        (Circular)                                               #topology
        \s*
        ([^.]+)                                                  #unknown
        \s*
        ($NPB::Parse::Format::Genbank::id_date)                             #date
        /xo) {

	$self->test_args($line, $1, $2, $3, $4, $5, $6, $7);

        (
         $self->{'id'},
         $self->{'length'},
         $self->{'units'},
         $self->{'molecule'},
         $self->{'topology'},
         $self->{'unknown'},
         $self->{'date'}

        ) = ($1, $2, $3, $4, $5, $6, $7);

    } else {
        $self->warn("unmatched: $line");
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n", 'id',       $self->{'id'};
    printf "$x%20s -> %s\n", 'length',   $self->{'length'};
    printf "$x%20s -> %s\n", 'units',    $self->{'units'};
    printf "$x%20s -> %s\n", 'molecule', $self->{'molecule'};
    printf "$x%20s -> %s\n", 'topology', $self->{'topology'};
    printf "$x%20s -> %s\n", 'unknown',  $self->{'unknown'};
    printf "$x%20s -> %s\n", 'date',     $self->{'date'};
}


###########################################################################
package NPB::Parse::Format::Genbank::ACCESSION;

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
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::Genbank::skip);

    $self->{'accessions'} = [];

    while (defined ($line = $text->next_line)) {
	$line =~ tr/;\n/ /;    #remove newlines and semicolon separators
	push @{$self->{'accessions'}}, split(" ", $line);
    }

    if (@{$self->{'accessions'}} < 1) {
	$self->warn("unmatched: $line");
    }

    $self->{'ac'} = $self->{'accessions'}->[0];
    
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n",   'ac',         $self->{'ac'};
    printf "$x%20s -> [%s]\n", 'accessions', join(',', @{$self->{'accessions'}});
}


###########################################################################
package NPB::Parse::Format::Genbank::FEATURES;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::FT);

#all definitions inherited


###########################################################################
package NPB::Parse::Format::Genbank::DEFINITION;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::die(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::Genbank::KEYWORDS;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::die(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::Genbank::SOURCE;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::die(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::Genbank::REFERENCE;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::die(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::Genbank::COMMENT;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::die(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::Genbank::BASECOUNT;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::die(shift, "new() not implemented");
}


###########################################################################
1;
