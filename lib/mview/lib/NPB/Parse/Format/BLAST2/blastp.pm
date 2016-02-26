# -*- perl -*-
# Copyright (C) 1996-2006 Nigel P. Brown
# $Id: blastp.pm,v 1.9 2005/12/12 20:42:48 brown Exp $

###########################################################################
package NPB::Parse::Format::BLAST2::blastp;

use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::HEADER);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::SEARCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::RANK);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::MATCH);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::MATCH::ALN);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::WARNING;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::WARNING);


###########################################################################
package NPB::Parse::Format::BLAST2::blastp::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
