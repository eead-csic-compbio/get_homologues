# -*- perl -*-
# Copyright (C) 1996-2013 Nigel P. Brown
# $Id: fasta.pm,v 1.10 2013/09/09 21:31:06 npb Exp $

###########################################################################
package NPB::Parse::Format::GCG_FASTA2::fasta;

use NPB::Parse::Format::GCG_FASTA2;
use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::GCG_FASTA2);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::GCG_FASTA2::fasta::HEADER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::GCG_FASTA2::HEADER);


###########################################################################
package NPB::Parse::Format::GCG_FASTA2::fasta::RANK;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::GCG_FASTA2::RANK);


###########################################################################
package NPB::Parse::Format::GCG_FASTA2::fasta::TRAILER;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::GCG_FASTA2::TRAILER);


###########################################################################
package NPB::Parse::Format::GCG_FASTA2::fasta::MATCH;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::GCG_FASTA2::MATCH);


###########################################################################
package NPB::Parse::Format::GCG_FASTA2::fasta::MATCH::SUM;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::GCG_FASTA2::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::GCG_FASTA2::fasta::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::GCG_FASTA2::MATCH::ALN);


###########################################################################
1;
