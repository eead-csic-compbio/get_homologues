# -*- perl -*-
# Copyright (C) 2011-2015 Nigel P. Brown
# $Id: fastx.pm,v 1.3 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::FASTA3::fastx;

use NPB::Parse::Format::FASTA3;
use NPB::Parse::Format::FASTA3::fasta;
use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::HEADER);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::RANK);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::TRAILER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::TRAILER);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::MATCH);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::fasta::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::FASTA3::fastx::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3::MATCH::ALN);

# fast[xy]  dna x pro

sub query_orient { $_[0]->get_summary->{'orient'} }
sub sbjct_orient { '+' }

sub query_base { return 3 }
sub sbjct_base { return 1 }


###########################################################################
1;
