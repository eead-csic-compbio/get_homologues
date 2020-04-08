# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST1::tblastn;

use Bio::Parse::Format::BLAST1::blastx;
use Bio::Parse::Format::BLAST1::tblastx;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::blastx);


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::blastx::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::tblastx::RANK);


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::blastx::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::blastx::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::MATCH::ALN;

use vars qw(@ISA);
use Bio::Util::Regexp;

@ISA = qw(Bio::Parse::Format::BLAST1::MATCH::ALN);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #Score line
    $line = $scan->read_line;

    if ($line =~ /^\s*
        Score\s*=\s*
        ($RX_Uint)                           #score
        \s+
        \(($RX_Ureal)\s+bits\),              #bits
        \s+
        Expect\s*=\s*
        ($RX_Ureal),                         #expectation
        \s+
        (?:Sum\sP\((\d+)\)|P)\s*=\s*         #number of frags
        ($RX_Ureal)                          #p-value
        /xo) {

        $self->test_args(\$line, $1, $2, $3, $5);

        (
         $self->{'score'},
         $self->{'bits'},
         $self->{'expect'},
         $self->{'n'},                       #substitute 1 unless $4
         $self->{'p'},
        ) = ($1, $2, $3, defined $4?$4:1, $5);
    }
    else {
        $self->warn("expecting 'Score' line: $line");
    }

    #Identities line
    $line = $scan->read_line;

    if ($line =~ /^\s*
        Identities\s*=\s*
        (\d+\/\d+)                           #identities fraction
        \s+
        \((\d+)%\),                          #identities percentage
        \s+
        Positives\s*=\s*
        (\d+\/\d+)                           #positives fraction
        \s+
        \((\d+)%\)                           #positives percentage
        ,\s+Frame\s*=\s*
        ([+-])                               #frame sign
        (\d+)                                #frame number
        /xo) {

        $self->test_args(\$line, $1, $2, $3, $4, $5, $6);

        (
         $self->{'id_fraction'},
         $self->{'id_percent'},
         $self->{'pos_fraction'},
         $self->{'pos_percent'},
         $self->{'sbjct_frame'},
        ) = ($1, $2, $3, $4, $5 . $6);

        #record sbjct orientation in MATCH list
        push @{$self->get_parent(1)->{'orient'}->{$self->{'sbjct_frame'}}}, $self;

    } else {
        $self->warn("expecting 'Identities' line: $line");
    }

    $self->parse_alignment($scan);

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = $self->SUPER::dump_data($indent);
    $s .= sprintf "$x%20s -> %s\n",  'score',        $self->{'score'};
    $s .= sprintf "$x%20s -> %s\n",  'bits',         $self->{'bits'};
    $s .= sprintf "$x%20s -> %s\n",  'expect',       $self->{'expect'};
    $s .= sprintf "$x%20s -> %s\n",  'p',            $self->{'p'};
    $s .= sprintf "$x%20s -> %s\n",  'n',            $self->{'n'};
    $s .= sprintf "$x%20s -> %s\n",  'id_fraction',  $self->{'id_fraction'};
    $s .= sprintf "$x%20s -> %s\n",  'id_percent',   $self->{'id_percent'};
    $s .= sprintf "$x%20s -> %s\n",  'pos_fraction', $self->{'pos_fraction'};
    $s .= sprintf "$x%20s -> %s\n",  'pos_percent',  $self->{'pos_percent'};
    $s .= sprintf "$x%20s -> %s\n",  'sbjct_frame',  $self->{'sbjct_frame'};
    return $s;
}


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::WARNING;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::WARNING);


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::HISTOGRAM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::HISTOGRAM);


###########################################################################
package Bio::Parse::Format::BLAST1::tblastn::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
