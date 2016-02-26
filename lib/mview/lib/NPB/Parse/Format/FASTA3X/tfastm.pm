# -*- perl -*-
# Copyright (C) 2015 Nigel P. Brown

###########################################################################
package NPB::Parse::Format::FASTA3X::tfastm;

use NPB::Parse::Format::FASTA3X;
use NPB::Parse::Format::FASTA3X::fastm;
use NPB::Parse::Format::FASTA3X::tfasta;

use strict;
use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3X);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package NPB::Parse::Format::FASTA3X::tfastm::HEADER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3X::fasta::HEADER);


###########################################################################
package NPB::Parse::Format::FASTA3X::tfastm::RANK;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3X::fastm::RANK);


###########################################################################
package NPB::Parse::Format::FASTA3X::tfastm::TRAILER;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA::TRAILER);


###########################################################################
package NPB::Parse::Format::FASTA3X::tfastm::MATCH;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA::MATCH);


###########################################################################
package NPB::Parse::Format::FASTA3X::tfastm::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(NPB::Parse::Format::FASTA3X::tfasta::MATCH::SUM);


###########################################################################
package NPB::Parse::Format::FASTA3X::tfastm::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(NPB::Parse::Format::FASTA3X::MATCH::ALN);

# tfast[mfs]  pro x dna

sub query_orient { '+' }
sub sbjct_orient { $_[0]->get_summary->{'orient'} }

sub query_base { return 1 }
sub sbjct_base { return 3 }

#Although this a S-W alignment, the sequence range numbers look more
#like those of a global/global alignment: query sequence range is complete;
#sbjct depends on it
sub get_start_stop {
    my ($self, $tgt, $orient, $info, $base) = @_;

    my ($start_info, $stop_info, $fs1, $fs2) = @$info;

    my ($r1) = @{$start_info};
    my ($r2) = @{$stop_info};

    warn "$tgt(o): $orient,$base start/stop: $r1,$r2\n\n"  if $NPB::Parse::Format::FASTA::DEBUG;
    return ($r1, $r2);
}


###########################################################################
1;
