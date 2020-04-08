# Copyright (C) 1996-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# FASTA 1
#
#   fasta, tfastx
#
###########################################################################
use Bio::MView::Build::Format::FASTA;

use strict;


###########################################################################
###########################################################################
package Bio::MView::Build::Row::FASTA1;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 1,   1,    'initn',         'initn',      '5N',      ''  ],
    [ 2,   2,    'init1',         'init1',      '5N',      ''  ],
    [ 3,   3,    'opt',           'opt',        '5N',      ''  ],
    [ 4,   4,    'query_orient',  'qy',         '2S',      '?' ],
    [ 5,   5,    'sbjct_orient',  'ht',         '2S',      '?' ],
    ]
}


###########################################################################
package Bio::MView::Build::Row::FASTA1::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA1);


###########################################################################
package Bio::MView::Build::Row::FASTA1::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA1);


###########################################################################
###########################################################################
package Bio::MView::Build::Format::FASTA1;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA);


###########################################################################
package Bio::MView::Build::Format::FASTA1::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA1);

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

    #warn "\nfasta1::parse: entering\n";

    #all strands done?
    if (! defined $self->{scheduler}->next) {
        $self->{'entry'}->free_parsers();
        #warn "fasta1::parse: exiting at end\n";
        return;
    }

    my $ranking = $self->{'entry'}->parse(qw(RANK));

    #no hits
    if (! defined $ranking) {
        #warn "fasta1::parse: exiting on empty parse\n";
        return;
    }

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
    my $class = "Bio::MView::Build::Row::FASTA1::$rtype";

    $coll->insert((new $class(
                       '',
                       $query,
                       '',
                       {
                           'initn'        => '',
                           'init1'        => '',
                           'opt'          => '',
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

        my $key1 = $coll->key($rank, $hit->{'id'}, $hit->{'initn'});

        $coll->insert((new $class(
                           $rank,
                           $hit->{'id'},
                           $hit->{'desc'},
                           {
                               'initn'        => $hit->{'initn'},
                               'init1'        => $hit->{'init1'},
                               'opt'          => $hit->{'opt'},
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

        my $key1 = $coll->key($rank, $sum->{'id'}, $sum->{'initn'});

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

            #for gapped alignments
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
        $coll->item($key1)->set_val('sbjct_orient', $sorient);
    }

    #warn "fasta1::parse: returning with data\n";

    return $coll->list;
}


###########################################################################
package Bio::MView::Build::Format::FASTA1::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA1);

sub parse {
    my $self = shift;

    #warn "\nfasta1::parse: entering\n";

    #all strands done?
    if (! defined $self->{scheduler}->next) {
        $self->{'entry'}->free_parsers();
        #warn "fasta1::parse: exiting at end\n";
        return;
    }

    my $ranking = $self->{'entry'}->parse(qw(RANK));

    #no hits
    if (! defined $ranking) {
        #warn "fasta1::parse: exiting on empty parse\n";
        return;
    }

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
    my $class = "Bio::MView::Build::Row::FASTA1::$rtype";

    $coll->insert((new $class(
                       '',
                       $query,
                       '',
                       {
                           'query_orient' => $self->strand,
                       }
                   )));

    #extract hits and identifiers from the ranking
    my $rank = 0; foreach my $hit (@{$ranking->{'hit'}}) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $hit->{'id'});

        #warn "KEEP: ($rank,$hit->{'id'})\n";

        my $key1 = $coll->key($rank, $hit->{'id'}, $hit->{'initn'});

        $coll->insert((new $class(
                           $rank,
                           $hit->{'id'},
                           $hit->{'desc'},
                           {
                               'initn'        => $hit->{'initn'},
                               'init1'        => $hit->{'init1'},
                               'opt'          => $hit->{'opt'},
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

        my $key1 = $coll->key($rank, $sum->{'id'}, $sum->{'initn'});

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
        $coll->item($key1)->set_val('sbjct_orient', $sorient);
    }

    #warn "fasta1::parse: returning with data\n";

    return $coll->list;
}


###########################################################################
1;
