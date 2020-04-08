# Copyright (C) 1996-2006 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST2::blastp;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::RANK);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::ALN);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::WARNING;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::WARNING);


###########################################################################
package Bio::Parse::Format::BLAST2::blastp::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
