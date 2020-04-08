# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST2::blastx;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2);


###########################################################################
package Bio::Parse::Format::BLAST2::blastx::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST2::blastx::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH);


###########################################################################
package Bio::Parse::Format::BLAST2::blastx::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::RANK);


###########################################################################
package Bio::Parse::Format::BLAST2::blastx::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST2::blastx::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST2::blastx::SEARCH::MATCH::ALN;

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

    #blastx 2.0.5 has no Frame line, so set a useful default
    $self->{'query_frame'} = $self->{'query_orient'};

    #query frame
    if ($line =~ /^\s*Frame\s*=\s+(\S+)\s*$/) {
        $self->{'query_frame'} = $1;
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
    $s .= sprintf "$x%20s -> %s\n",  'query_frame',  $self->{'query_frame'};
    return $s;
}

###########################################################################
package Bio::Parse::Format::BLAST2::blastx::WARNING;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::WARNING);


###########################################################################
package Bio::Parse::Format::BLAST2::blastx::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
