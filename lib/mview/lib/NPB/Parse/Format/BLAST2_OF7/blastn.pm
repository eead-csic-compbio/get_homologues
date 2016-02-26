# -*- perl -*-
# Copyright (C) 2015 Nigel P. Brown
# $Id: blastn.pm $

###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastn;

use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastn::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::HEADER);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastn::SEARCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastn::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::RANK);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastn::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::MATCH);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastn::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::BLAST2_OF7::blastn::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN);

sub new {
    my $type = shift;
    my $self = new NPB::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN(@_);

    #record paired orientations in MATCH list
    push @{$self->get_parent(1)->{'orient'}->{
				 $self->{'query_orient'} .
				 $self->{'sbjct_orient'}
				}}, $self;
    bless $self, $type;
}

###########################################################################
1;
