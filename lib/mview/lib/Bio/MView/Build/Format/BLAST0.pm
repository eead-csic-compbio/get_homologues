# Copyright (C) 2015-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# NCBI BLAST 1.4, WashU BLAST 2.0
#
#   blastp, blastn, blastx, tblastn, tblastx
#
# NCBI BLAST 2.0, PSI-BLAST 2.0, BLAST+
#
#   blastp, blastn, blastx, tblastn, tblastx
#
###########################################################################
package Bio::MView::Build::Format::BLAST0;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Format::BLAST;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST);

#override
sub begin_parse {}   #sets up parser objects
sub end_parse {}     #frees parser objects
sub parse_record {}  #parse data for a row record
sub record_type {}   #returns type of row to instantiate
sub skip_hsp {}      #score/signifiance filter

#override/extend
sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s  if $quiet;
    my $mode = $PAR->get('hsp');
    if ($mode eq 'all') {
        $s .= "HSP processing: all\n";
    } elsif ($mode eq 'discrete') {
        $s .= "HSP processing: discrete\n";
    } else {
        $s .= "HSP processing: ranked\n";
    }
    $s;
}

sub hi_score {
    my ($self, $score1, $aln) = @_;
    my $score2 = $aln->{$self->{'attr_score'}};
    return $score2 if $score2 > $score1;
    $score1;
}

sub lo_sig {
    my ($self, $sig1, $aln) = @_;
    my $sig2 = $aln->{$self->{'attr_sig'}};
    return $sig2 if $sig2 < $sig1 or $sig1 < 0;
    $sig1;
}

sub get_scores {
    my ($self, $list) = @_;

    return unless @$list;

    my ($score, $qorient, $sorient) = (0, '?', '?');

    foreach my $aln (@$list) {
        $score = $aln->{$self->{'attr_score'}}  if
            $aln->{$self->{'attr_score'}} > $score;
    }

    foreach my $aln (@$list) {
        $qorient = $aln->{'query_orient'}, next  if $qorient eq '?';
        if ($aln->{'query_orient'} ne $qorient) {
            warn "get_scores: mixed up query orientations\n";
        }
    }

    foreach my $aln (@$list) {
        $sorient = $aln->{'sbjct_orient'}, next  if $sorient eq '?';
        if ($aln->{'sbjct_orient'} ne $sorient) {
            warn "get_scores: mixed up sbjct orientations\n";
        }
    }

    my $n = scalar @$list;
    my $sig = $list->[0]->{$self->{'attr_sig'}};

    ($n, $score, $sig, $qorient, $sorient);
}

sub parse {
    my $self = shift;

    #warn "\nblast0::parse: entering\n";

    if (! defined $self->{scheduler}->next) {
        $self->end_parse;
        #warn "blast0::parse: exiting at end\n";
        return;
    }

    my ($search, $header, $ranking) = $self->begin_parse;

    #no hits
    if (! ($search and $ranking)) {
        #warn "blast0::parse: exiting on empty parse\n";
        return;
    }

    my $coll = new Bio::MView::Build::Search::Collector($self);

    $coll->insert($self->record(undef, $header, undef, undef));

    #extract hits and identifiers from the ranking
    my $rank = 0; foreach my $hit (@{$ranking->{'hit'}}) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $hit->{'id'});

        #warn "KEEP: ($rank,$hit->{'id'})\n";

        my $key1 = $coll->key($hit->{'id'});

        $coll->insert($self->record($rank, undef, undef, $hit), $key1);
    }

    my $mode = $PAR->get('hsp');

    if ($mode eq 'all') {
        $self->parse_all_hits($coll, $search);

    } elsif ($mode eq 'ranked') {
        $self->parse_ranked_hits($coll, $search);

    } elsif ($mode eq 'discrete') {
        $self->parse_discrete_hits($coll, $search);
    }

    #warn "blast0::parse: returning with data\n";

    return $coll->list;
}

sub parse_ranked_hits {
    my ($self, $coll, $search) = @_;

    #pull out each hit
    foreach my $match ($search->parse(qw(MATCH))) {

        #first the summary
        my $sum = $match->parse(qw(SUM));

        my $key1 = $coll->key($sum->{'id'});

        #ignore hit?
        next  unless $coll->has($key1);

        #$self->report_ranking_data($match, $coll, $key1, $self->strand), next;

        my $raln = $self->get_ranked_hsps($match, $coll, $key1, $self->strand);

        #nothing matched
        next  unless @$raln;

        foreach my $aln (@$raln) {
            #apply score/significance filter
            next  if $self->skip_hsp($aln);

            #for gapped alignments: keep sbjct sequence insertions?
            my $uncut = ($PAR->get('keepinserts') ? $aln->{'sbjct'} : '');
            $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'});

            $coll->add_frags(
                $key1, $aln->{'query_start'}, $aln->{'query_stop'}, [
                    $aln->{'query'},
                    $aln->{'query_start'},
                    $aln->{'query_stop'},
                    $aln->{$self->{'attr_score'}},
                ], [
                    $aln->{'sbjct'},
                    $aln->{'sbjct_start'},
                    $aln->{'sbjct_stop'},
                    $aln->{$self->{'attr_score'}},
                    \$uncut,
                ]);
        }
        #override row data
        $coll->item($key1)->{'desc'} = $sum->{'desc'};

        my ($N, $score, $sig, $qorient, $sorient) = $self->get_scores($raln);
        $coll->item($key1)->set_val('n', $N);
        $coll->item($key1)->set_val($self->{'attr_score'}, $score);
        $coll->item($key1)->set_val($self->{'attr_sig'}, $sig);
        $coll->item($key1)->set_val('query_orient', $qorient);
        $coll->item($key1)->set_val('sbjct_orient', $sorient);
    }
    $self;
}

sub parse_discrete_hits {
    my ($self, $coll, $search) = @_;

    #pull out each hit
    foreach my $match ($search->parse(qw(MATCH))) {

        #first the summary
        my $sum = $match->parse(qw(SUM));

        my $key1 = $coll->key($sum->{'id'});

        #ignore hit?
        next  unless $coll->has($key1);

        foreach my $aln ($match->parse(qw(ALN))) {

            #ignore other query strand orientation
            next  unless $aln->{'query_orient'} eq $self->strand;

            my $key2 = $coll->key($match->{'index'}, $aln->{'index'});

            #apply row filter with new row numbers
            next  if $self->skip_row($match->{'index'}, $key2, $sum->{'id'});
            next  if $self->skip_hsp($aln);

            #make new row
            if (! $coll->has($key2)) {
                $coll->insert($self->record($key2, undef, $sum, $aln), $key2);
            }

            #for gapped alignments: keep sbjct sequence insertions?
            my $uncut = ($PAR->get('keepinserts') ? $aln->{'sbjct'} : '');
            $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'});

            $coll->add_frags(
                $key2, $aln->{'query_start'}, $aln->{'query_stop'}, [
                    $aln->{'query'},
                    $aln->{'query_start'},
                    $aln->{'query_stop'},
                    $aln->{$self->{'attr_score'}},
                ], [
                    $aln->{'sbjct'},
                    $aln->{'sbjct_start'},
                    $aln->{'sbjct_stop'},
                    $aln->{$self->{'attr_score'}},
                    \$uncut,
                ]);

            #override row data
            $coll->item($key2)->set_val('query_orient', $aln->{'query_frame'})
                if exists $aln->{'query_frame'};
            $coll->item($key2)->set_val('sbjct_orient', $aln->{'sbjct_frame'})
                if exists $aln->{'sbjct_frame'};

            #override N: for discrete output this must always be 1
            $coll->item($key2)->set_val('n', 1);
        }
    }
    $self;
}

sub parse_blastpx_all_hits {
    my ($self, $coll, $search) = @_;

    #pull out each hit
    foreach my $match ($search->parse(qw(MATCH))) {

        #first the summary
        my $sum = $match->parse(qw(SUM));

        my $key1 = $coll->key($sum->{'id'});

        #ignore hit?
        next  unless $coll->has($key1);

        my ($n, $score, $sig) = (0, 0, -1);

        foreach my $aln ($match->parse(qw(ALN))) {

            #ignore other query strand orientation
            next  unless $aln->{'query_orient'} eq $self->strand;

            #apply score/significance filter
            next  if $self->skip_hsp($aln);

            #accumulate row data
            $score = $self->hi_score($score, $aln);
            $sig   = $self->lo_sig($sig, $aln);
            $n++;

            #for gapped alignments: keep sbjct sequence insertions?
            my $uncut = ($PAR->get('keepinserts') ? $aln->{'sbjct'} : '');
            $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'});

            $coll->add_frags(
                $key1, $aln->{'query_start'}, $aln->{'query_stop'}, [
                    $aln->{'query'},
                    $aln->{'query_start'},
                    $aln->{'query_stop'},
                    $aln->{$self->{'attr_score'}},
                ], [
                    $aln->{'sbjct'},
                    $aln->{'sbjct_start'},
                    $aln->{'sbjct_stop'},
                    $aln->{$self->{'attr_score'}},
                    \$uncut,
                ]);
        }
        #override row data
        $coll->item($key1)->{'desc'} = $sum->{'desc'};
        $coll->item($key1)->set_val('n', $n);
        $coll->item($key1)->set_val($self->{'attr_score'}, $score);
        $coll->item($key1)->set_val($self->{'attr_sig'}, $sig);
    }
    $self;
}


sub parse_tblastnx_all_hits {
    my ($self, $coll, $search) = @_;

    #pull out each hit
    foreach my $match ($search->parse(qw(MATCH))) {

        #first the summary
        my $sum = $match->parse(qw(SUM));

        my $key1 = $coll->key($sum->{'id'});

        #ignore hit?
        next  unless $coll->has($key1);

        my ($n1,$n2, $score1,$score2, $sig1,$sig2) = (0,0,  0,0, -1,-1);

        foreach my $aln ($match->parse(qw(ALN))) {

            #ignore other query strand orientation
            next  unless $aln->{'query_orient'} eq $self->strand;

            #apply score/significance filter
            next  if $self->skip_hsp($aln);

            #accumulate row data
            my $rank   = $coll->key($match->{'index'}, $aln->{'index'});

            if ($aln->{'sbjct_orient'} eq '+') {
                $score1 = $self->hi_score($score1, $aln);
                $sig1   = $self->lo_sig($sig1, $aln);
                $n1++;
            } else {
                $score2 = $self->hi_score($score2, $aln);
                $sig2   = $self->lo_sig($sig2, $aln);
                $n2++;
            }

            my $key2 = $coll->key($sum->{'id'}, $aln->{'sbjct_orient'});

            if (! $coll->has($key2)) {
                $coll->insert($self->record($rank, undef, $sum, $aln), $key2);
            }

            #for gapped alignments: keep sbjct sequence insertions?
            my $uncut = ($PAR->get('keepinserts') ? $aln->{'sbjct'} : '');
            $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'});

            $coll->add_frags(
                $key2, $aln->{'query_start'}, $aln->{'query_stop'}, [
                    $aln->{'query'},
                    $aln->{'query_start'},
                    $aln->{'query_stop'},
                    $aln->{$self->{'attr_score'}},
                ], [
                    $aln->{'sbjct'},
                    $aln->{'sbjct_start'},
                    $aln->{'sbjct_stop'},
                    $aln->{$self->{'attr_score'}},
                    \$uncut,
                ]);
        }
        #override row data
        $coll->item($key1)->{'desc'} = $sum->{'desc'};

        my $keyp = $coll->key($key1, '+');
        if ($coll->has($keyp)) {
            $coll->item($keyp)->set_val('n', $n1);
            $coll->item($keyp)->set_val($self->{'attr_score'}, $score1);
            $coll->item($keyp)->set_val($self->{'attr_sig'}, $sig1);
        }

        my $keym = $coll->key($key1, '-');
        if ($coll->has($keym)) {
            $coll->item($keym)->set_val('n', $n2);
            $coll->item($keym)->set_val($self->{'attr_score'}, $score2);
            $coll->item($keym)->set_val($self->{'attr_sig'}, $sig2);
        }
    }
    $self;
}


###########################################################################
package Bio::MView::Build::Format::BLAST0::blastp;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST0);

sub record {
    my $self = shift;
    my $type = $self->record_type;
    my @values = $self->parse_record(@_);

    $values[3]->{'query_orient'} = '+';
    $values[3]->{'sbjct_orient'} = '+';

    new $type(@values);
}

sub parse_all_hits { my $self = shift; $self->parse_blastpx_all_hits(@_) }


###########################################################################
package Bio::MView::Build::Format::BLAST0::blastn;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST0);

sub record {
    my $self = shift;
    my $type = $self->record_type;
    my @values = $self->parse_record(@_);

    my $qorient = $self->strand;
    my $sorient = '?';

    if (defined $_[3]) {
        $sorient = $_[3]->{'sbjct_orient'}  if exists $_[3]->{'sbjct_orient'};
    } else {
        $sorient = ''; #must be the query row: has no sbjct orientation
    }

    $values[3]->{'query_orient'} = $qorient;
    $values[3]->{'sbjct_orient'} = $sorient;

    new $type(@values);
}

sub parse_all_hits { my $self = shift; $self->parse_tblastnx_all_hits(@_) }


###########################################################################
package Bio::MView::Build::Format::BLAST0::blastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST0);

sub record {
    my $self = shift;
    my $type = $self->record_type;
    my @values = $self->parse_record(@_);

    my $qorient = $self->strand;

    $values[3]->{'query_orient'} = $qorient;
    $values[3]->{'sbjct_orient'} = '+';

    new $type(@values);
}

sub parse_all_hits { my $self = shift; $self->parse_blastpx_all_hits(@_) }


###########################################################################
package Bio::MView::Build::Format::BLAST0::tblastn;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST0);

sub record {
    my $self = shift;
    my $type = $self->record_type;
    my @values = $self->parse_record(@_);

    my $sorient = '?';

    if (defined $_[3]) {
        $sorient = $_[3]->{'sbjct_orient'}  if exists $_[3]->{'sbjct_orient'};
    } else {
        $sorient = ''; #must be the query row: has no sbjct orientation
    }

    $values[3]->{'query_orient'} = '+';
    $values[3]->{'sbjct_orient'} = $sorient;

    new $type(@values);
}

sub parse_all_hits { my $self = shift; $self->parse_tblastnx_all_hits(@_) }


###########################################################################
package Bio::MView::Build::Format::BLAST0::tblastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST0);

sub record {
    my $self = shift;
    my $type = $self->record_type;
    my @values = $self->parse_record(@_);

    my $qorient = $self->strand;
    my $sorient = '?';

    if (defined $_[3]) {
        $qorient = $_[3]->{'query_orient'}  if exists $_[3]->{'query_orient'};
        $sorient = $_[3]->{'sbjct_orient'}  if exists $_[3]->{'sbjct_orient'};
    } else {
        $sorient = ''; #must be the query row: has no sbjct orientation
    }

    $values[3]->{'query_orient'} = $qorient;
    $values[3]->{'sbjct_orient'} = $sorient;

    new $type(@values);
}

sub parse_all_hits { my $self = shift; $self->parse_tblastnx_all_hits(@_) }


###########################################################################
