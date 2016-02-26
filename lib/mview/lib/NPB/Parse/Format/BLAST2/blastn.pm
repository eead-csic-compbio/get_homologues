# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: blastn.pm,v 1.11 2015/06/14 17:09:04 npb Exp $

###########################################################################
package NPB::Parse::Format::BLAST2::blastn;

use strict;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2);


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::HEADER);


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::SEARCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH);


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::RANK);


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::MATCH);


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST2::SEARCH::MATCH::ALN);

sub new {
    my $type = shift;
    my ($parent) = @_;
    my $self = new NPB::Parse::Format::BLAST2::SEARCH::MATCH::ALN(@_);
    bless $self, $type;

    #record paired orientations in MATCH list
    push @{$parent->{'orient'}->{
				 $self->{'query_orient'} .
				 $self->{'sbjct_orient'}
				}}, $self;
    $self;
}


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::WARNING;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::WARNING);


###########################################################################
package NPB::Parse::Format::BLAST2::blastn::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
