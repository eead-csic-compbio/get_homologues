# Copyright (C) 2015-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastx;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastx::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastx::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastx::SEARCH::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::RANK);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastx::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastx::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::blastx::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN);

my $MAP_ALN = { 'qframe' => 'query_frame' };

sub new {
    my $type = shift;
    my $self = new Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN(@_);
    my $text = new Bio::Parse::Scanner($self);

    #BLAST2 -outfmt 7
    $self->{'query_frame'} = '';

    $self->extract_fields($MAP_ALN, $text->read_line(1));

    if ($self->{'query_frame'} eq '') {
        $self->die("blast column specifier 'qframe' is needed");
    }

    #prepend sign to forward frame
    $self->{'query_frame'} = "+$self->{'query_frame'}"
        if $self->{'query_frame'} =~ /^\d+/;

    #record paired orientations in MATCH list
    push @{$self->get_parent(1)->{'orient'}->{
               $self->{'query_orient'} . $self->{'sbjct_orient'}
           }}, $self;

    bless $self, $type;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = $self->SUPER::dump_data($indent);
    $s .= sprintf "$x%20s -> %s\n",  'query_frame',  $self->{'query_frame'};
    return $s;
}


###########################################################################
1;
