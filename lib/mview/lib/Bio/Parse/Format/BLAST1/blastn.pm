# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::BLAST1::blastn;

use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1);

my $NULL      = '^\s*$';
my $RANK_NONE = '^\s*\*\*\* NONE';
my $RANK_CUT  = 60;


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::RANK;

use Bio::Parse::Strings qw(strip_trailing_space clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #column headers
    $self->{'header'} = $scan->read_lines(4);

    #ranked search hits
    $self->{'hit'}    = [];

    while (defined ($line = $scan->read_line(1))) {

        #blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

        #GCG annotation: ignore
        next    if $line =~ /$Bio::Parse::Format::BLAST::GCG_JUNK/o;

        #empty ranking: done
        last    if $line =~ /$RANK_NONE/o;

        #excise variable length description and append it
        my $tmp = substr($line, 0, $RANK_CUT);
        if ($tmp =~ /^\s*([^\s]+)(.*)/o) {
            $line = $1 . substr($line, $RANK_CUT) . $2;
        }

        if ($line =~ /\s*
            ([^\s]+)                          #id
            \s+
            ($RX_Uint)                        #score
            \s+
            ($RX_Ureal)                       #p-value
            \s+
            ($RX_Uint)                        #n fragments
            \s*
            \!?                               #GCG junk
            \s*
            (.*)                              #summary
            /xo) {

            $self->test_args(\$line, $1, $2, $3, $4); #ignore $5

            push @{$self->{'hit'}},
            {
             'id'      => clean_identifier($1),
             'score'   => $2,
             'p'       => $3,
             'n'       => $4,
             'summary' => strip_trailing_space($5),
            };

            next;
        }

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST1::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::MATCH::ALN;

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
        ($RX_Uint)                      #score
        \s+
        \(($RX_Ureal)\s+bits\),         #bits
        \s+
        Expect\s*=\s*
        ($RX_Ureal),                    #expectation
        \s+
        (?:Sum\sP\((\d+)\)|P)\s*=\s*    #number of frags
        ($RX_Ureal)                     #p-value
        /xo) {

        $self->test_args(\$line, $1, $2, $3, $5);

        (
         $self->{'score'},
         $self->{'bits'},
         $self->{'expect'},
         $self->{'n'},                  #substitute 1 unless $4
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
        (\d+\/\d+)                      #identities fraction
        \s+
        \((\d+)%\),                     #identities percentage
        \s+
        Positives\s*=\s*
        (\d+\/\d+)                      #positives fraction
        \s+
        \((\d+)%\)                      #positives percentage
        ,\s+Strand\s+=\s+
        (\S+)\s+\/\s+(\S+)              #strand orientations
        /xo) {

        $self->test_args(\$line, $1, $2, $3, $4, $5, $6);

        (
         $self->{'id_fraction'},
         $self->{'id_percent'},
         $self->{'pos_fraction'},
         $self->{'pos_percent'},
         $self->{'query_orient'},
         $self->{'sbjct_orient'},
        ) = ($1, $2, $3, $4, $5 eq 'Plus'?'+':'-', $6 eq 'Plus'?'+':'-');

        #record query orientation in MATCH list
        push @{$self->get_parent(1)->{'orient'}->{$self->{'query_orient'}}}, $self;

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
    return $s;
}


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::WARNING;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::WARNING);


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::HISTOGRAM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::HISTOGRAM);


###########################################################################
package Bio::Parse::Format::BLAST1::blastn::PARAMETERS;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::PARAMETERS);


###########################################################################
1;
