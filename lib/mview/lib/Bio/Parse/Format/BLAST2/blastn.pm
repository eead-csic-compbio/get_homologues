# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST2::blastn;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2);


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH);


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::RANK);


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::ALN);

sub new {
    my $type = shift;
    my $self = new Bio::Parse::Format::BLAST2::SEARCH::MATCH::ALN(@_);
    bless $self, $type;

    #record paired orientations in MATCH list
    push @{$self->get_parent(1)->{'orient'}->{
               $self->{'query_orient'} . $self->{'sbjct_orient'}
           }}, $self;

    $self;
}


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::WARNING;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::WARNING);


###########################################################################
package Bio::Parse::Format::BLAST2::blastn::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
