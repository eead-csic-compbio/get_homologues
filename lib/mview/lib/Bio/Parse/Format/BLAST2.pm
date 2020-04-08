# Copyright (C) 1996-2020 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Base classes for NCBI BLAST2 family.
#
# Handles: BLAST2 and BLAST+
#
# BLAST (NCBI version 2) iterated searching uses 3 main record types:
#
#   HEADER        the header text
#   SEARCH {*}    passes of the search engine
#     RANK        the list of ordered high scoring hits
#     MATCH {*}   the set of alignments for a given hit
#       SUM         the summary lines for each hit
#       ALN {*}     each aligned fragment: score + alignment
#   PARAMETERS    the trailer
#
###########################################################################
package Bio::Parse::Format::BLAST2;

use Bio::Parse::Format::BLAST;
use Bio::Util::Regexp;

use strict;

use vars qw(@ISA

            @VERSIONS

            $ENTRY_START
            $ENTRY_END

            $WARNING_START
            $WARNING_END

            $WARNINGS_START
            $WARNINGS_END

            $PARAMETERS_START
            $PARAMETERS_END

            $MATCH_START
            $MATCH_END

            $RANK_START
            $RANK_MATCH
            $RANK_END
            $RANK_NONE

            $HEADER_START
            $HEADER_END

            $SEARCH_START
            $SEARCH_END

            $SCORE_START
            $SCORE_END
           );

@ISA   = qw(Bio::Parse::Format::BLAST);

@VERSIONS = (
             '2' => [
                     'BLASTP',
                     'BLASTN',
                     'BLASTX',
                     'TBLASTN',
                     'TBLASTX',
                     'PSIBLAST',
                     'PHIBLASTP',
                    ],
            );

my $NULL = '^\s*$';

$ENTRY_START = '^(?:'
    . 'BLASTP'
    . '|'
    . 'BLASTN'
    . '|'
    . 'BLASTX'
    . '|'
    . 'TBLASTN'
    . '|'
    . 'TBLASTX'
    . '|'
    . 'PSIBLAST'
    . '|'
    . 'PHIBLASTP'
    . ')';
$ENTRY_END        = '^  WARNINGS\s+ISSUED:';  #blast1 behavior, but blast2?

$SEARCH_START     = "^(?:Searching|Results from round|\\s+High\\s+E|\\s+Score\\s+E)";
$SEARCH_END       = "^(?:$SEARCH_START|  Database|CPU time:|Lambda +)";

$WARNING_START    = '^WARNING';
$WARNING_END      = "(?:done|$NULL)";  #'done' keyword in phi-blast SEARCH block

$WARNINGS_START   = '^WARNINGS\s+ISSUED:';
$WARNINGS_END     = $NULL;

$PARAMETERS_START = "^(?:  Database|CPU time:|Lambda +)";
$PARAMETERS_END   = "^(?:$WARNINGS_START|$ENTRY_START|$SEARCH_START)";

$MATCH_START      = '^>';
$MATCH_END        = "^(?:$MATCH_START|$WARNING_START|$SEARCH_END)";

$RANK_START       = "^(?:\\s+High\\s+E|\\s+Score\\s+E)";
$RANK_MATCH       = "^(?:$NULL|\\s+High|\\s+Score|\\s*Sequences|\\s*CONVERGED| *No hits found|[^>].*$RX_Uint\\s+$RX_Ureal|$Bio::Parse::Format::BLAST::GCG_JUNK)";
$RANK_END         = "(?:$MATCH_START|$SEARCH_END)";  #ignore warnings
$RANK_NONE        = '^\s*\*\*\* NONE';

$HEADER_START     = $ENTRY_START;
$HEADER_END       = "^(?:$SEARCH_START|$WARNING_START|$RANK_START|$MATCH_START|$PARAMETERS_START)";

$SCORE_START      = '^\s{0,1}Score';  #exactly 0 or 1 leading spaces
$SCORE_END        = "^(?:$SCORE_START|$MATCH_END)";


sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

        #Header lines
        if ($line =~ /$HEADER_START/o) {
            $scan->scan_until($HEADER_END);
                $self->push_record('HEADER',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #Search lines
        if ($line =~ /$SEARCH_START/o) {
            if ($line =~ /^Searching/o) {  #BLAST2
                $scan->scan_until("^(?:Searching|$PARAMETERS_START)");
                $self->push_record('SEARCH',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
                next;
            }
            if ($line =~ /^Results from round/o) {  #BLAST+
                $scan->scan_until("^(?:Results from round|$PARAMETERS_START)");
                $self->push_record('SEARCH',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
                next;
            }
            #keywords 'Searching' or 'Results' stripped by web server
            $scan->scan_until($SEARCH_END);
                $self->push_record('SEARCH',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #Parameter lines
        if ($line =~ /$PARAMETERS_START/o) {
            $scan->scan_until($PARAMETERS_END);
                $self->push_record('PARAMETERS',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #WARNINGS ISSUED line: ignore
        next    if $line =~ /$WARNINGS_START/o;

        #default
        $self->warn("unknown field: $line");
    }

    $self;#->examine;
}

#BLAST2 or PSI-BLAST write  'e10' instead of '1e10', so fix this.
sub fix_expect {
    my $e = shift;
    $e =~ s/^e/1e/;
    $e;
}


###########################################################################
package Bio::Parse::Format::BLAST2::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::HEADER);


###########################################################################
package Bio::Parse::Format::BLAST2::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line(1))) {
        #warn "[$line]\n";

        #query line: update HEADER if needed (for psiblast at least)
        if ($line =~ /^Query=\s+(\S+)?\s*[,;]?\s*(.*)/o) {
            #no test - either field may be missing
            my $query   = $1  if defined $1;
            my $summary = $2  if defined $2;

            #strip leading unix path stuff
            $query =~ s/.*\/([^\/]+)$/$1/;

            my $header = $self->get_parent(1)->get_record('HEADER');
            $header->{'query'}   = $query    if $header->{'query'} eq '';
            $header->{'summary'} = $summary  if $header->{'summary'} eq '';
        }

        #Rank lines
        if ($line =~ /$Bio::Parse::Format::BLAST2::RANK_START/o) {
            $scan->scan_until($Bio::Parse::Format::BLAST2::RANK_END);
                $self->push_record('RANK',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #Hit lines
        if ($line =~ /$Bio::Parse::Format::BLAST2::MATCH_START/o) {
            $scan->scan_until($Bio::Parse::Format::BLAST2::MATCH_END);
                $self->push_record('MATCH',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #WARNING lines
        if ($line =~ /$Bio::Parse::Format::BLAST2::WARNING_START/o) {
            $scan->scan_until($Bio::Parse::Format::BLAST2::WARNING_END);
                $self->push_record('WARNING',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #default: ignore any other text in SEARCH block
    }

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::BLAST2::SEARCH::RANK;

use Bio::Parse::Strings qw(strip_trailing_space clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #column headers
    $self->{'header'} = $scan->read_until_inclusive('Value(?:\s*N)?\s*$');

    #ranked search hits
    $self->{'hit'}    = [];

    while (defined ($line = $scan->read_line(1))) {

        next    if $line =~ /^Sequences used in model and found again:/o;
        next    if $line =~ /^Sequences not found previously or not previously below threshold:/o;
        next    if $line =~ /^CONVERGED!/o;
        next    if $line =~ /^Significant /o;    #PHI-BLAST

        #blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

        #GCG annotation: ignore
        next    if $line =~ /$Bio::Parse::Format::BLAST::GCG_JUNK/o;

        #NCBI web server HTML text-only output: ignore
        next    if $line =~ /^\s*-+\s*$/o;        #dashed line
        next    if $line =~ /^\s+Alignments/o;    #text

        #empty ranking: done
        last    if $line =~ /$Bio::Parse::Format::BLAST2::RANK_NONE/o;

        my $tmp = {};

        #examine suffix
        if ($line =~ /
            \s+
            ($RX_Ureal)                #bits
            \s+
            ($RX_Ureal)                #E value
            (?:\s+($RX_Uint))?         #N (ungapped blast only)
            \s*$
            /xo) {

            $self->test_args(\$line, $1, $2);

            $tmp->{'bits'}   = $1;
            $tmp->{'expect'} = Bio::Parse::Format::BLAST2::fix_expect($2);
            $tmp->{'n'}      = (defined $3 ? $3 : 0);

            #examine prefix
            if ($` =~ /^\s*
                (\S+)                  #id
                \s*
                \!?                    #GCG junk
                \s*
                (.*)?                  #summary
                /xo) {

                $self->test_args(\$line, $1);    #ignore $2

                $tmp->{'id'} = clean_identifier($1);
                $tmp->{'summary'} = strip_trailing_space($2);

            } else {
                $self->warn("unknown field: $line");
            }

            push @{$self->{'hit'}}, $tmp;

            next;
        }

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::BLAST2::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::MATCH);


###########################################################################
package Bio::Parse::Format::BLAST2::SEARCH::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::BLAST2::SEARCH::MATCH::ALN;

use vars qw(@ISA);
use Bio::Util::Regexp;

@ISA = qw(Bio::Parse::Format::BLAST::MATCH::ALN);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #Score line
    $line = $scan->read_line;

    if ($line =~ /^\s*
        Score\s*=\s*
        ($RX_Ureal)                 #bits
        \s+bits\s+
        \(($RX_Ureal)\),            #score
        \s+
        Expect(?:\((\d+)\))?\s*=\s*
        ($RX_Ureal)                 #expectation
        /xo) {

        $self->test_args(\$line, $1, $2, $4);

        (
         $self->{'bits'},
         $self->{'score'},
         $self->{'n'},              #substitute 1 unless $3
         $self->{'expect'},
        ) = ($1, $2, defined $3?$3:1, Bio::Parse::Format::BLAST2::fix_expect($4));

    } else {
        $self->warn("expecting 'Score' line: $line");
    }

    #Identities line
    $line = $scan->read_line;

    if ($line =~ /^\s*
        Identities\s*=\s*
        (\d+\/\d+)                  #identities fraction
        \s+
        \((\d+)%\)                  #identities percentage
        (?:                         #not present in (some?) BLASTN 2.0.9
        ,\s+
        Positives\s*=\s*
        (\d+\/\d+)                  #positives fraction
        \s+
        \((\d+)%\)                  #positives percentage
        (?:,\s*                     #not always present
         Gaps\s*=\s*
         (\d+\/\d+)                 #gaps fraction
         \s+
         \((\d+)%\)                 #gaps percentage
        )?
        )?
        /xo) {

        $self->test_args(\$line, $1, $2);

        (
         $self->{'id_fraction'},
         $self->{'id_percent'},
         $self->{'pos_fraction'},
         $self->{'pos_percent'},
         $self->{'gap_fraction'},
         $self->{'gap_percent'},
        ) = ($1, $2, defined $3?$3:'', defined $4?$4:0, defined $5?$5:'', defined $6?$6:0);

    } else {
        $self->warn("expecting 'Identities' line: $line");
    }

    #optional (Strand=|Frame=) line: handled in subclasses
    $line = $scan->read_line;

    $self->parse_alignment($scan);

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = $self->SUPER::dump_data($indent);
    $s .= sprintf "$x%20s -> %s\n", 'bits',         $self->{'bits'};
    $s .= sprintf "$x%20s -> %s\n", 'score',        $self->{'score'};
    $s .= sprintf "$x%20s -> %s\n", 'n',            $self->{'n'};
    $s .= sprintf "$x%20s -> %s\n", 'expect',       $self->{'expect'};
    $s .= sprintf "$x%20s -> %s\n", 'id_fraction',  $self->{'id_fraction'};
    $s .= sprintf "$x%20s -> %s\n", 'id_percent',   $self->{'id_percent'};
    $s .= sprintf "$x%20s -> %s\n", 'pos_fraction', $self->{'pos_fraction'};
    $s .= sprintf "$x%20s -> %s\n", 'pos_percent',  $self->{'pos_percent'};
    $s .= sprintf "$x%20s -> %s\n", 'gap_fraction', $self->{'gap_fraction'};
    $s .= sprintf "$x%20s -> %s\n", 'gap_percent',  $self->{'gap_percent'};
    return $s;
}


###########################################################################
1;
