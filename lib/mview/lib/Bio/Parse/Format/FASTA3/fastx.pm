# Copyright (C) 2011-2015 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA3::fastx;

use Bio::Parse::Format::FASTA3;
use Bio::Parse::Format::FASTA3::fasta;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA3::fastx::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3::fasta::HEADER);


###########################################################################
package Bio::Parse::Format::FASTA3::fastx::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3::fasta::RANK);


###########################################################################
package Bio::Parse::Format::FASTA3::fastx::TRAILER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3::fasta::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA3::fastx::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3::fasta::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA3::fastx::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3::fasta::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::FASTA3::fastx::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3::MATCH::ALN);

# fast[xy]  dna x pro

sub query_orient { $_[0]->get_summary->{'orient'} }
sub sbjct_orient { '+' }

sub query_base { return 3 }
sub sbjct_base { return 1 }


###########################################################################
1;
