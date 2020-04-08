# Copyright (C) 1996-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# FASTA 2
#
#   fasta, tfastx (tested)
#   fastx, fasty, tfasta, tfasty, tfastxy (not tested: may work)
#
###########################################################################
use Bio::MView::Build::Format::FASTA;

use strict;


###########################################################################
###########################################################################
package Bio::MView::Build::Row::FASTA2;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 1,   1,    'initn',         'initn',      '5N',      ''  ],
    [ 2,   2,    'init1',         'init1',      '5N',      ''  ],
    [ 3,   3,    'opt',           'opt',        '5N',      ''  ],
    [ 4,   4,    'zscore',        'z-sc',       '7N',      ''  ],
    [ 5,   5,    'expect',        'E-value',    '9N',      ''  ],
    [ 6,   6,    'query_orient',  'qy',         '2S',      '?' ],
    [ 7,   7,    'sbjct_orient',  'ht',         '2S',      '?' ],
    ]
}


###########################################################################
package Bio::MView::Build::Row::FASTA2::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA2);


###########################################################################
package Bio::MView::Build::Row::FASTA2::tfasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA2::fasta);


###########################################################################
package Bio::MView::Build::Row::FASTA2::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA2::fasta);


###########################################################################
package Bio::MView::Build::Row::FASTA2::tfasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA2::fasta);


###########################################################################
package Bio::MView::Build::Row::FASTA2::tfastxy;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA2::fasta);


###########################################################################
###########################################################################
package Bio::MView::Build::Format::FASTA2;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA);


###########################################################################
package Bio::MView::Build::Format::FASTA2::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2);

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Query orientation: " . $self->strand . "\n";
    $s;
}

sub parse {
    my $self = shift;

    #warn "\nfasta2::parse: entering\n";

    #all strands done?
    if (! defined $self->{scheduler}->next) {
        $self->{'entry'}->free_parsers();
        #warn "fasta2::parse: exiting at end\n";
        return;
    }

    #identify the query
    my $header = $self->{'entry'}->parse(qw(HEADER));

    my $ranking = $self->{'entry'}->parse(qw(RANK));

    #no hits
    if (! defined $ranking) {
        #warn "fasta2::parse: exiting on empty parse\n";
        return;
    }

    my $query = 'Query';
    if ($header->{'query'} ne '') {
        $query = $header->{'query'};
    } elsif ($header->{'queryfile'} =~ m,.*/([^\.]+)\.,) {
        $query = $1;
    }

    my $coll = new Bio::MView::Build::Search::Collector($self);

    my $rtype = $1  if ref($self) =~ /::([^:]+)$/;
    my $class = "Bio::MView::Build::Row::FASTA2::$rtype";

    $coll->insert((new $class(
                       '',
                       $query,
                       '',
                       {
                           'initn'        => '',
                           'init1'        => '',
                           'opt'          => '',
                           'zscore'       => '',
                           'expect'       => '',
                           'query_orient' => $self->strand,
                           'sbjct_orient' => '',
                       }
                   )));

    #extract hits and identifiers from the ranking
    my $rank = 0; foreach my $hit (@{$ranking->{'hit'}}) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $hit->{'id'});

        #warn "KEEP: ($rank,$hit->{'id'})\n";

        my $key1 = $coll->key($rank, $hit->{'id'}, $hit->{'expect'});

        $coll->insert((new $class(
                           $rank,
                           $hit->{'id'},
                           $hit->{'desc'},
                           {
                               'initn'        => $hit->{'initn'},
                               'init1'        => $hit->{'init1'},
                               'opt'          => $hit->{'opt'},
                               'zscore'       => $hit->{'zscore'},
                               'expect'       => $hit->{'expect'},
                               'query_orient' => $self->strand,
                               'sbjct_orient' => '',
                           }
                       )), $key1
            );
    }

    #pull out each hit
    $rank = 0; foreach my $match ($self->{'entry'}->parse(qw(MATCH))) {

        $rank++;

        #first the summary
        my $sum = $match->parse(qw(SUM));

        my $key1 = $coll->key($rank, $sum->{'id'}, $sum->{'expect'});

        next  unless $coll->has($key1);

        #warn "USE: [$key1]\n";

        my $aset = $self->get_ranked_frags($match, $coll, $key1, $self->strand);

        #nothing matched
        next  unless @$aset;

        foreach my $aln (@$aset) {
            #apply score/significance filter
            next  if $self->skip_frag($sum->{'opt'});

            #warn "SEE: [$key1]\n";

            #$aln->dump;

            #for FASTA gapped alignments
            $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'},
                                    $aln->{'query_leader'},
                                    $aln->{'query_trailer'});

            $coll->add_frags(
                $key1, $aln->{'query_start'}, $aln->{'query_stop'}, [
                    $aln->{'query'},
                    $aln->{'query_start'},
                    $aln->{'query_stop'},
                ], [
                    $aln->{'sbjct'},
                    $aln->{'sbjct_start'},
                    $aln->{'sbjct_stop'},
                ]);
        }
        #override row data
        $coll->item($key1)->{'desc'} = $sum->{'desc'}  if $sum->{'desc'};

        my ($N, $sig, $qorient, $sorient) = $self->get_scores($aset, $sum);
        #$coll->item($key1)->set_val('n', $N);        #not used yet
        $coll->item($key1)->set_val('init1', $sum->{'init1'});
        $coll->item($key1)->set_val('initn', $sum->{'initn'});
        #$coll->item($key1)->set_val('expect', $sig); #not used yet
        $coll->item($key1)->set_val('sbjct_orient', $sorient);
    }

    #warn "fasta2::parse: returning with data\n";

    return $coll->list;
}


###########################################################################
package Bio::MView::Build::Format::FASTA2::fastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2::fasta);


###########################################################################
package Bio::MView::Build::Format::FASTA2::fasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2::fasta);


###########################################################################
package Bio::MView::Build::Format::FASTA2::tfasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2);

#note 2015-07-04: for future refactoring: this differs unusually from
#fasta::parse because the sbjct orientation has to be derived from the SUM
#instead of the ALN object, i.e., ALN orientations derived from the sequence
#numbering are wrong.
sub parse {
    my $self = shift;

    #all strands done?
    return  unless defined $self->{scheduler}->next;

    my $ranking = $self->{'entry'}->parse(qw(RANK));

    #fasta run with no hits
    return []  unless defined $ranking;

    #identify the query
    my $header = $self->{'entry'}->parse(qw(HEADER));

    my $query = 'Query';
    if ($header->{'query'} ne '') {
        $query = $header->{'query'};
    } elsif ($header->{'queryfile'} =~ m,.*/([^\.]+)\.,) {
        $query = $1;
    }

    my $coll = new Bio::MView::Build::Search::Collector($self);

    my $rtype = $1  if ref($self) =~ /::([^:]+)$/;
    my $class = "Bio::MView::Build::Row::FASTA2::$rtype";

    $coll->insert((new $class(
                       '',
                       $query,
                       '',
                       {
                           'initn'        => '',
                           'init1'        => '',
                           'opt'          => '',
                           'zscore'       => '',
                           'expect'       => '',
                           'query_orient' => $self->strand,
                           'sbjct_orient' => '',
                       }
                   )));

    #extract hits and identifiers from the ranking
    my $rank = 0; foreach my $hit (@{$ranking->{'hit'}}) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $hit->{'id'});

        #warn "KEEP: ($rank,$hit->{'id'})\n";

        my $key1 = $coll->key($rank, $hit->{'id'}, $hit->{'expect'});

        $coll->insert((new $class(
                           $rank,
                           $hit->{'id'},
                           $hit->{'desc'},
                           {
                               'initn'        => $hit->{'initn'},
                               'init1'        => $hit->{'init1'},
                               'opt'          => $hit->{'opt'},
                               'zscore'       => $hit->{'zscore'},
                               'expect'       => $hit->{'expect'},
                               'query_orient' => $self->strand,
                               'sbjct_orient' => $hit->{'orient'},
                           }
                       )), $key1
            );
    }

    #pull out each hit
    $rank = 0; foreach my $match ($self->{'entry'}->parse(qw(MATCH))) {

        $rank++;

        #first the summary
        my $sum = $match->parse(qw(SUM));

        my $key1 = $coll->key($rank, $sum->{'id'}, $sum->{'expect'});

        next  unless $coll->has($key1);

        #warn "USE: [$key1]\n";

        my $aset = $self->get_ranked_frags($match, $coll, $key1, $self->strand);

        #nothing matched
        next  unless @$aset;

        foreach my $aln (@$aset) {
            #apply score/significance filter
            next  if $self->skip_frag($sum->{'opt'});

            #warn "SEE: [$key1]\n";

            #$aln->dump;

            #for FASTA gapped alignments
            $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'},
                                    $aln->{'query_leader'},
                                    $aln->{'query_trailer'});

            $coll->add_frags(
                $key1, $aln->{'query_start'}, $aln->{'query_stop'}, [
                    $aln->{'query'},
                    $aln->{'query_start'},
                    $aln->{'query_stop'},
                ], [
                    $aln->{'sbjct'},
                    $aln->{'sbjct_start'},
                    $aln->{'sbjct_stop'},
                ]);

            #override row data (non-standard)
            $coll->item($key1)->set_val('sbjct_orient', $sum->{'orient'});
        }
        #override row data
        $coll->item($key1)->{'desc'} = $sum->{'desc'}  if $sum->{'desc'};

        my ($N, $sig, $qorient, $sorient) = $self->get_scores($aset, $sum);
        #$coll->item($key1)->set_val('n', $N);        #not used yet
        $coll->item($key1)->set_val('init1', $sum->{'init1'});
        $coll->item($key1)->set_val('initn', $sum->{'initn'});
        #$coll->item($key1)->set_val('expect', $sig); #not used yet
        #$coll->item($key1)->set_val('sbjct_orient', $sorient); NO: see above
    }

    #free objects
    $self->{'entry'}->free_parsers(qw(HEADER RANK MATCH));

    return $coll->list;
}


###########################################################################
package Bio::MView::Build::Format::FASTA2::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2::tfasta);


###########################################################################
package Bio::MView::Build::Format::FASTA2::tfasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2::tfasta);


###########################################################################
package Bio::MView::Build::Format::FASTA2::tfastxy;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2::tfasta);


###########################################################################
###########################################################################
package Bio::MView::Build::Format::FASTA2::ssearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2);

sub parse {
    warn "FASTA2 SSEARCH processing - not implemented\n";
    return undef;
}


###########################################################################
package Bio::MView::Build::Format::FASTA2::align;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2);

sub parse {
    warn "FASTA2 ALIGN processing - not implemented\n";
    return undef;
}


###########################################################################
package Bio::MView::Build::Format::FASTA2::lalign;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA2);

sub parse {
    warn "FASTA2 LALIGN processing - not implemented\n";
    return undef;
}


###########################################################################
1;
