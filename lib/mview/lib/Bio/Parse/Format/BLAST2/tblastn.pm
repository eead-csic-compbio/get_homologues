# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST2::tblastn;

use Bio::Parse::Format::BLAST2::blastx;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::blastx);


###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::blastx::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::blastx::SEARCH);


###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::blastx::SEARCH::RANK);


###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::blastx::SEARCH::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::blastx::SEARCH::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::ALN);

sub new {
    my $type = shift;
    my $self = new Bio::Parse::Format::BLAST2::SEARCH::MATCH::ALN(@_);
    bless $self, $type;

    my $scan = new Bio::Parse::Scanner($self);
    my $line;

    $line = $scan->read_line(1);
    $line = $scan->read_line(1);
    $line = $scan->read_line(1);

    #warn "[$line]\n";

    #tblastn 2.0.5 has no Frame line, so set a useful default
    $self->{'sbjct_frame'} = $self->{'sbjct_orient'};

    #sbjct frame
    if ($line =~ /^\s*Frame\s*=\s+(\S+)\s*$/o) {
        $self->{'sbjct_frame'} = $1;
    }

    #record paired orientations in MATCH list
    push @{$self->get_parent(1)->{'orient'}->{
               $self->{'query_orient'} . $self->{'sbjct_orient'}
           }}, $self;

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = $self->SUPER::dump_data($indent);
    $s .= sprintf "$x%20s -> %s\n",  'sbjct_frame',  $self->{'sbjct_frame'};
    return $s;
}

###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::WARNING;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::WARNING);


###########################################################################
package Bio::Parse::Format::BLAST2::tblastn::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
