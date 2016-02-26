# -*- perl -*-
# Copyright (C) 2015 Nigel P. Brown
# $Id: blastp.pm $

###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastp;

use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastp::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::HEADER);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastp::SEARCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastp::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::RANK);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastp::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::MATCH);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastp::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastp::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN);


###########################################################################
1;
