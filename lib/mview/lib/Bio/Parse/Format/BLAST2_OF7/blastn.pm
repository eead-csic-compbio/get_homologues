# Copyright (C) 2015-2017 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastn;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastn::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastn::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastn::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::RANK);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastn::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastn::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastn::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN);

sub new {
    my $type = shift;
    my $self = new Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN(@_);

    #record paired orientations in MATCH list
    push @{$self->get_parent(1)->{'orient'}->{
               $self->{'query_orient'} . $self->{'sbjct_orient'}
           }}, $self;

    bless $self, $type;
}

###########################################################################
1;
