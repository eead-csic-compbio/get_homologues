# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA2_GCG::fasta;

use Bio::Parse::Format::FASTA2_GCG;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA2_GCG);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::fasta::HEADER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2_GCG::HEADER);


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::fasta::RANK;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2_GCG::RANK);


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::fasta::TRAILER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2_GCG::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::fasta::MATCH;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2_GCG::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::fasta::MATCH::SUM;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2_GCG::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::fasta::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2_GCG::MATCH::ALN);


###########################################################################
1;
