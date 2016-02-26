# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: EMBL.pm,v 1.52 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::EMBL;

use NPB::Parse::FT;
use vars qw(@ISA);
use strict;

@ISA = qw(NPB::Parse::Record);


#Exported constants
#skip these many chars at front of each raw record
$NPB::Parse::Format::EMBL::skip          = 5;

#special EMBL record fields
$NPB::Parse::Format::EMBL::id_entryname  = '[A-Z]\w{1,8}';
$NPB::Parse::Format::EMBL::ac_accession  = $NPB::Parse::FT::ac_accession;


#EMBL record types
my $EMBL_ID            = '^ID';    
my $EMBL_AC            = '^AC';    
my $EMBL_NI            = '^NI';    
my $EMBL_DT            = '^DT';    
my $EMBL_DE            = '^DE';    
my $EMBL_KW            = '^KW';    
my $EMBL_OS            = '^OS';    
my $EMBL_OG            = '^OG';    
my $EMBL_RN            = '^RN';    
my $EMBL_DR            = '^DR';    
my $EMBL_CC            = '^CC';    
my $EMBL_FH            = '^FH';    
my $EMBL_FT            = '^FT';    
my $EMBL_SQ            = '^SQ';    
my $EMBL_XX            = '^XX';        # end of simple or compound record
my $EMBL_End           = '^\/\/';      # end of entry: two slashes
my $EMBL_Null          = '^\S*\s*$';#' # blank line or empty record {FH,KW}


#Consume one entry-worth of input on stream $fh associated with $file and
#return a new EMBL instance.
sub get_entry {
    my ($parent) = @_;
    my ($line, $offset, $bytes) = ('', -1, 0);

    my $fh   = $parent->{'fh'};
    my $text = $parent->{'text'};

    while (defined ($line = <$fh>)) {
	
	#ID line: start of entry
	if ($line =~ /$EMBL_ID/o and $offset < 0) {
	    $offset = $fh->tell - length($line);
	    next;
	}

	#// line: end of entry
	if ($line =~ /$EMBL_End/o) {
	    last;
	}
    }
    return 0   if $offset < 0;

    $bytes = $fh->tell - $offset;

    new NPB::Parse::Format::EMBL(undef, $text, $offset, $bytes);
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

	#ID line
	if ($line =~ /$EMBL_ID/o) {
	    $record = $text->scan_while($EMBL_ID, 'ID');
	    $ob = new NPB::Parse::Format::EMBL::ID($self, \$record);
	    $self->push_indices($ob->{'id'});
	    next;
	}

	#AC line
	if ($line =~ /$EMBL_AC/o) {
	    $record = $text->scan_while($EMBL_AC, 'AC');
	    $ob = new NPB::Parse::Format::EMBL::AC($self, \$record);
	    $self->set_relative_key($ob->{'ac'});
	    ##$self->push_indices(@{$ob->{'accessions'}});
	    next;
	}

	#NI line
	if ($line =~ /$EMBL_NI/o) {
	    $record = $text->scan_while($EMBL_AC, 'NI');
	    $ob = new NPB::Parse::Format::EMBL::NI($self, \$record);
	    ##$self->push_indices(@{$ob->{'accessions'}});
	    next;
	}

	#DT line
	if ($line =~ /$EMBL_DT/o) {
	    $text->scan_while($EMBL_DT, 'DT');
	    next;			       	      
	}				       	      
	
	#DE line			       	      
	if ($line =~ /$EMBL_DE/o) {       	      
	    $text->scan_while($EMBL_DE, 'DE');
	    next;			       	      
	}				       	      
	
	#KW line			       	      
	if ($line =~ /$EMBL_KW/o) {       	      
	    $text->scan_while($EMBL_KW, 'KW');
	    next;			       	      
	}				       	      
	
	#OS line (+ lines: OC)
	if ($line =~ /$EMBL_OS/o) {       	      
	    $text->scan_until($EMBL_XX, 'OS');
	    next;
	}
	
	#OG line			       	      
	if ($line =~ /$EMBL_OG/o) {       	      
	    $text->scan_while($EMBL_OG, 'OG');
	    next;
	}
	
	#RN line (+lines: RP, RX, RA, RT, RL, RC)
	if ($line =~ /$EMBL_RN/o) {		     
	    $text->scan_until($EMBL_XX, 'RN');
	    next;					     
	}						     
	
	#DR line			             	     
	if ($line =~ /$EMBL_DR/o) {             	     
	    $text->scan_while($EMBL_DR, 'DR');
	    next;			             	     
	}				             	     
	
	#CC line			             	     
	if ($line =~ /$EMBL_CC/o) {             	     
	    $text->scan_while($EMBL_CC, 'CC');
	    next;			             	     
	}				             	     
	
	#FH line: ignore
	if ($line =~ /$EMBL_FH/o) {
	    next;
	} 
	
	#FT line			             	     
	if ($line =~ /$EMBL_FT/o) {             	     
	    $text->scan_while($EMBL_FT, 'FT');
	    next;			             	     
	}				             	     
	
	#SQ line 
	if ($line =~ /$EMBL_SQ/o) {             	     
	    $text->scan_until($EMBL_End, 'SQ');
	    next;
	}
	
	#XX line: ignore
	if ($line =~ /$EMBL_XX/o) {
	    next;
	} 
	
	#blank line or empty record: ignore
	if ($line =~ /$EMBL_Null/o) {
	    next;
	}

	#terminal line: ignore
	if ($line =~ /$EMBL_End/o) {
	    next;
	}

	#default
	$self->warn("unknown field: ", substr($line, 0, $NPB::Parse::Format::EMBL::skip));
    }
    $self;
}


###########################################################################
package NPB::Parse::Format::EMBL::ID;

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
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::EMBL::skip);

    $line = $text->next_line;
    
    if ($line =~ /^\s*
        ($NPB::Parse::Format::EMBL::id_entryname)                     #id
        \s*
        (standard|unreviewed|preliminary|unannotated|backbone)   #dataclass
        \s*;\s* 
        (DNA|RNA|xxx|circular\ DNA|circular\ RNA|circular\ xxx)  #molecule
        \s*;\s*
        ([A-Z]{3})                                               #division
        \s*;\s*
        ($RX_Uint)                                               #length
        \s*
        ([^.]+)\.                                                #units
        /xo) {

        $self->test_args($line, $1, $2, $3, $4, $5, $6);

        (
         $self->{'id'},
         $self->{'dataclass'},
         $self->{'molecule'},
         $self->{'division'},
         $self->{'length'},
         $self->{'units'}

        ) = ($1, $2, $3, $4, $5, $6);

    } else {
        $self->warn("unmatched: $line");
    }
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n", 'id',        $self->{'id'};
    printf "$x%20s -> %s\n", 'dataclass', $self->{'dataclass'};
    printf "$x%20s -> %s\n", 'molecule',  $self->{'molecule'};
    printf "$x%20s -> %s\n", 'division',  $self->{'division'};
    printf "$x%20s -> %s\n", 'length',    $self->{'length'};
    printf "$x%20s -> %s\n", 'units',     $self->{'units'};
}


###########################################################################
package NPB::Parse::Format::EMBL::AC;

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
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::EMBL::skip);

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
package NPB::Parse::Format::EMBL::NI;

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
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::EMBL::skip);

    $line = $text->next_line;

    if ($line =~ /^\s*
	(\w+)                                                    #g14412
	/xo) {

	$self->test_args($line, $1);

	$self->{'ni'} = $1;
	
    }

    if ($self->{'ni'} eq '') {
	$self->warn("unmatched: $line");
    }

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> %s\n", 'ni', $self->{'ni'};
}


###########################################################################
package NPB::Parse::Format::EMBL::DT;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::warn(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::EMBL::DE;

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
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::EMBL::skip);

    $line = $text->scan_lines(0);

    $self->{'description'} = NPB::Parse::Record::strip_english_newlines($line);

    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> \"%s\"\n", 'description', $self->{'description'};
}


###########################################################################
package NPB::Parse::Format::EMBL::FT;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::FT);

#all definitions inherited


###########################################################################
package NPB::Parse::Format::EMBL::KW;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::warn(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::EMBL::OS;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::warn(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::EMBL::OG;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::warn(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::EMBL::RN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::warn(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::EMBL::DR;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::warn(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::EMBL::CC;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Record);

sub new {
    NPB::Message::warn(shift, "new() not implemented");
}


###########################################################################
package NPB::Parse::Format::EMBL::SQ;

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
    $text = new NPB::Parse::Record_Stream($self, $NPB::Parse::Format::EMBL::skip);

    $self->{'sequence'} = '';

    while (defined ($line = $text->next_line)) {

	#sequence header line
	if ($line =~ /^\s*
	    Sequence
	    \s+
	    (\d+)                                                #length
	    \s+BP;
	    \s+(\d+)\sA;                                         #A
	    \s+(\d+)\sC;                                         #C
	    \s+(\d+)\sG;                                         #G
	    \s+(\d+)\sT;                                         #T
	    \s+(\d+)\sother;                                     #other
	    \s*
	    /xo) {

	    $self->test_args($line, $1, $2, $3, $4, $5, $6);

	    (
	     $self->{'length'},
	     $self->{'A'},
	     $self->{'C'},
	     $self->{'G'},
	     $self->{'T'},
	     $self->{'other'}
	     
	    ) = ($1, $2, $3, $4, $5, $6);
	    
	    next;
	} 

	#sequence lines
	if ($line =~ /^\s*
	    (.*)                                                 #sequence
	    \s+
	    \d+
	    \s*
	    /xo) {

	    $self->test_args($line, $1);

	    $self->{'sequence'} .= $1;

	    next;
	}
	
    }

    $self->{'sequence'} =~ tr/ \n//d;    #strip space
    $self;
}

sub print_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    printf "$x%20s -> \"%s\"\n",  'sequence', $self->{'sequence'};
    printf "$x%20s -> %d\n",      'length',   $self->{'length'};
    printf "$x%20s -> %d\n",      'A',        $self->{'A'};
    printf "$x%20s -> %d\n",      'C',        $self->{'C'};
    printf "$x%20s -> %d\n",      'G',        $self->{'G'};
    printf "$x%20s -> %d\n",      'T',        $self->{'T'};
    printf "$x%20s -> %d\n",      'other',    $self->{'other'};
}


###########################################################################
1;
