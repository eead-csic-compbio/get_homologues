# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# generic BLAST material
#
###########################################################################
package Bio::MView::Build::Format::BLAST;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Search;
use Bio::Util::Regexp;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Search);

#the name of the underlying Bio::Parse::Stream parser
sub parser { 'BLAST' }

my $MISSING_QUERY_CHAR = 'X';  #interpolate this between query fragments

my %Known_Parameters =
    (
     #name        => [ format       default  ]

     #BLAST* display various HSP selections
     'hsp'         => [ '\S+',       'ranked' ],
     'keepinserts' => [ '\d+',              0 ],

     #BLAST* (version 1)
     'maxpval'     => [ $RX_Ureal,   undef    ],
     'minscore'    => [ $RX_Ureal,   undef    ],

     #BLAST* (version 2)
     'maxeval'     => [ $RX_Ureal,   undef    ],
     'minbits'     => [ $RX_Ureal,   undef    ],
     'cycle'       => [ [],          undef    ],

     #BLASTN (version 1, version 2); BLASTX (version 2)
     'strand'      => [ [],          undef    ],

     #BLASTX/TBLASTX (version 1)
    );

#standard attribute names: override
sub attr_score { 'unknown' }
sub attr_sig   { 'unknown' }

#scheduler type: override
sub scheduler { 'unknown' }

#our own constructor since this is the entry point for different subtypes
sub new {
    shift;  #discard type
    my $self = new Bio::MView::Build::Search(@_);
    my ($type, $p, $v, $file);

    #determine the real type from the underlying parser
    ($p, $v) = (lc $self->{'entry'}->{'format'}, $self->{'entry'}->{'version'});

    $type = "Bio::MView::Build::Format::BLAST$v";
    ($file = $type) =~ s/::/\//g;
    require "$file.pm";

    $type .= "::$p";
    bless $self, $type;

    $self->{'attr_score'} = $self->attr_score;
    $self->{'attr_sig'}   = $self->attr_sig;

    $self->initialise;

    $self;
}

#called by the constructor
sub initialise {
    my $self = shift;
    my $scheduler = $self->scheduler;
    #warn "initialise ($scheduler)\n";
    while (1) {
        $self->{scheduler} = new Bio::MView::Build::Scheduler,
        last if $scheduler eq 'none';

        $self->{scheduler} = new Bio::MView::Build::Scheduler([qw(+ -)]);
        last if $scheduler eq 'strand';

        if ($scheduler eq 'cycle') {
            my $last = $self->{'entry'}->count(qw(SEARCH));
            $self->{scheduler} = new Bio::MView::Build::Scheduler([1..$last]);
            last;
        }

        if ($scheduler eq 'cycle+strand') {
            my $last = $self->{'entry'}->count(qw(SEARCH));
            $self->{scheduler} =
                new Bio::MView::Build::Scheduler([1..$last], [qw(+ -)]);
            last;
        }

        die "initialise: unknown scheduler '$scheduler'";
    }
    return $self;
}

#called on each iteration
sub reset_child {
    my $self = shift;
    my $scheduler = $self->scheduler;
    #warn "reset_child ($scheduler)\n";

    my $strand = $PAR->get('strand');
    my $cycle  = $PAR->get('cycle');

    while (1) {
        last if $scheduler eq 'none';

        #(warn "strands: [@{$strand}]\n"),
        $self->{scheduler}->filter($strand),
        last if $scheduler eq 'strand';

        #(warn "cycles: [@{$cycle}]\n"),
        $self->{scheduler}->filter($cycle),
        last if $scheduler eq 'cycle';

        #(warn "cycles+strands: [@{$cycle}][@{$strand}]\n"),
        $self->{scheduler}->filter($cycle, $strand),
        last if $scheduler eq 'cycle+strand';

        die "reset_child: unknown scheduler '$scheduler'";
    }
    return $self;
}

#current cycle being processed
sub cycle {
    my $scheduler = $_[0]->scheduler;
    return $_[0]->{scheduler}->item       if $scheduler eq 'cycle';
    return ($_[0]->{scheduler}->item)[0]  if $scheduler eq 'cycle+strand';
    return 1;
}

#current strand being processed
sub strand {
    my $scheduler = $_[0]->scheduler;
    return $_[0]->{scheduler}->item       if $scheduler eq 'strand';
    return ($_[0]->{scheduler}->item)[1]  if $scheduler eq 'cycle+strand';
    return '+';
}

#override base class method to process query row differently
sub build_rows {
    my $self = shift;
    my ($lo, $hi) = (0,0);

    #compute alignment length from query sequence in row[0]
    if (! $PAR->get('keepinserts')) {
        ($lo, $hi) = $self->get_range($self->{'index2row'}->[0]);
    }
    #warn "BLAST::build_rows range[0] ($lo, $hi)\n";

    #query row contains missing query sequence, rather than gaps
    $self->{'index2row'}->[0]->assemble($lo, $hi, $MISSING_QUERY_CHAR);
    #assemble sparse sequence strings for all rows
    for (my $i=1; $i < @{$self->{'index2row'}}; $i++) {
        #warn "BLAST::build_rows range[$i] ($lo, $hi)\n";
        $self->{'index2row'}->[$i]->assemble($lo, $hi, $PAR->get('gap'));
    }
    $self;
}

#utility/testing function
sub report_ranking_data {
    my ($self, $match, $coll, $rkey, $qorient) = @_;
    return  unless $qorient eq '+';  #called multiply; only want one pass
    my $n     = $coll->item($rkey)->get_val('n');
    my $score = $coll->item($rkey)->get_val($self->{'attr_score'});
    my $sig   = $coll->item($rkey)->get_val($self->{'attr_sig'});
    my $state = -1;
    my ($asco, $asig);
    #look for a match (ranking==hit by score and sig) in either orientation
    foreach my $aln ($match->parse(qw(ALN))) {
        last  if $state > 1;
        $state = 0; #entered loop body at least once
        $asco = $aln->{$self->{'attr_score'}};
        $asig = $aln->{$self->{'attr_sig'}};
        #match conditions
        $state = 2, next  if $score == $asco and $sig == $asig;
        my $ascor = sprintf('%0.f', $asco);
        $state = 3, next  if $score == $ascor and $sig == $asig;
    }
    return  if $state < 0;
    warn "(@{[$self->cycle]})$rkey match: identity\n"  if $state == 2;
    warn "(@{[$self->cycle]})$rkey match: round(score=$score/$asco)\n"
        if $state == 3;
    return  if $state > 0;
    #no match, start reporting
    warn "(@{[$self->cycle]})$rkey <<<<<<<<<<<<<<<<\n";
    foreach my $aln ($match->parse(qw(ALN))) {
        my $asco = $aln->{$self->{'attr_score'}};
        my $asig = $aln->{$self->{'attr_sig'}};
        my $aqorient = $aln->{'query_orient'};
        warn "$aqorient score: $score $asco " . "sig: $sig $asig  " .
            "@{[$score == $asco ? '1' : '0']} @{[$sig == $asig ? '1' : '0']} | @{[$score == sprintf('%0.f', $asco) ? '1' : '0']}\n";
    }
}

#return a hash of HSP sets; each set is a list of ALN with a given query/sbjct
#ordering, number of fragments, and significance; these sets have multiple
#alternative keys: "qorient/sorient/n/sig", "qorient/sorient/n/score",
#"qorient/sorient/score" to handle fuzziness and errors in the input data.
sub get_hsp_groups {
    my ($self, $match) = @_;
    my $hash = {};
    foreach my $aln ($match->parse(qw(ALN))) {
        my $n       = $aln->{'n'};
        my $score   = $aln->{$self->{'attr_score'}};
        my $sig     = $aln->{$self->{'attr_sig'}};
        my $qorient = $aln->{'query_orient'};
        my $sorient = $aln->{'sbjct_orient'};
        my (%tmp, $key);

        #first ALN has precedence: we use this later for key hints
        $hash->{'.'} = $aln  if keys %$hash < 1;

        #key by significance: can exceed $n if more HSP have the same sig
        %tmp = ( $sig => 1 );               #raw significance value
        #BLAST1 and WU-BLAST
        if ($self->{'attr_sig'} eq 'p' and $sig !~ /e/i) {
            $tmp{sprintf("%.2f", $sig)}++;  #rounded sig (2dp)
        }
        #generic
        foreach my $sig (keys %tmp) {
            #full information
            $key = join("/", $qorient, $sorient, $n, $sig);
            push @{ $hash->{$key} }, $aln;

            #without N as this can be missing in the ranking or wrong
            $key = join("/", $qorient, $sorient, $sig);
            push @{ $hash->{$key} }, $aln;
        }

        #key by N and score, but we also need to consider rounding:
        #
        #round to nearest integer, or, if X.5, round both ways; this is
        #because the scores reported in the ranking and in the hit summary are
        #themselves rounded from an unknown underlying score, so simply
        #rounding the HSP score is not enough, viz.
        # underlying: 9.49 -> ranking (0dp) 9
        # underlying: 9.49 -> hit     (1dp) 9.5
        # round(9.5, 0) = 10 != 9
        %tmp = ( $score => 1 );             #raw score
        $tmp{sprintf("%.0f", $score)}++;    #rounded score (0dp)
        if ($score - int($score) == 0.5) {  #edge cases +/- 0.5
            $tmp{$score - 0.5}++;
            $tmp{$score + 0.5}++;
        }
        foreach my $score (keys %tmp) {
            #full information
            $key = join("/", $qorient, $sorient, $n, $score);
            push @{ $hash->{$key} }, $aln;

            #without N as this can be missing in the ranking or wrong
            $key = join("/", $qorient, $sorient, $score);
            push @{ $hash->{$key} }, $aln;
        }
    }
    #warn "SAVE: @{[sort keys %$hash]}\n";
    return $hash;
}

#lookup an HSP set in an HSP dictionary, trying various keys based on number
#of HSPs N, significance or score.
sub get_ranked_hsps_by_query {
    my ($self, $hash, $coll, $rkey, $qorient1) = @_;

    my $n        = $coll->item($rkey)->get_val('n');
    my $score    = $coll->item($rkey)->get_val($self->{'attr_score'});
    my $qorient2 = ($qorient1 eq '+' ? '-' : '+');
    my $sorient1  = $hash->{'.'}->{'sbjct_orient'}; #first ALN orientation hint
    my $sorient2 = ($sorient1 eq '+' ? '-' : '+');  #its opposite

    $n = 1  unless $n;  #no N in ranking

    my $D = 0;

    my $lookup_o_n_sig = sub {
        my ($qorient, $n, $sig) = @_;
        foreach my $sorient ($sorient1, $sorient2) {
            my $key = join("/", $qorient, $sorient, $n, $sig);
            return $hash->{$key}  if exists $hash->{$key};
        }
        return undef;
    };

    my $follow_sig = sub {
        my ($qorient, $key) = @_;

        #take the sig of the first matching HSP: expand that set
        my $n     = $hash->{$key}->[0]->{'n'};
        my $score = $hash->{$key}->[0]->{$self->{'attr_score'}};
        my $sig   = $hash->{$key}->[0]->{$self->{'attr_sig'}};
        warn "MAP: $qorient $n $score -> $sig\n"  if $D;

        my $match;

        #match in query orientation?
        $match = &$lookup_o_n_sig($qorient, $n, $sig);
        return $match  if defined $match;

        #match in opposite orientation? skip
        $match = &$lookup_o_n_sig($qorient2, $n, $sig);
        return []  if defined $match;

        return undef;
    };

    my $lookup_o_n_score = sub {
        my ($qorient, $n, $score) = @_;
        foreach my $sorient ($sorient1, $sorient2) {
            my $key = join("/", $qorient, $sorient, $n, $score);
            return &$follow_sig($qorient, $key)  if exists $hash->{$key};
        }
        return undef;
    };

    my $lookup_score = sub {
        my ($qorient, $score) = @_;
        foreach my $sorient ($sorient1, $sorient2) {
            my $key = join("/", $qorient, $sorient, $score);
            return &$follow_sig($qorient, $key)  if exists $hash->{$key};
        }
        return undef;
    };

    my $match;

    warn "KEYS($rkey,$qorient1): @{[sort keys %$hash]}\n"  if $D;

    #match (n, score) in query orientation?
    warn "TRY: $qorient1/$n/$score\n"  if $D;
    $match = &$lookup_o_n_score($qorient1, $n, $score);
    if (defined $match) {
        warn "MATCH (@{[scalar @$match]})\n"  if $D;
        return $match;
    }

    #match (n, score) in opposite orientation? skip
    warn "TRY: $qorient2/$n/$score (to skip)\n"  if $D;
    $match = &$lookup_o_n_score($qorient2, $n, $score);
    if (defined $match) {
        warn "SKIP (@{[scalar @$match]})\n"  if $D;
        return [];
    }

    #match (score) in query orientation?
    warn "TRY: $qorient1/$score\n"  if $D;
    $match = &$lookup_score($qorient1, $score);
    if (defined $match) {
        warn "MATCH (@{[scalar @$match]})\n"  if $D;
        return $match;
    }

    #match (score) in opposite orientation? skip
    warn "TRY: $qorient2/$score (to skip)\n"  if $D;
    $match = &$lookup_score($qorient2, $score);
    if (defined $match) {
        warn "SKIP (@{[scalar @$match]})\n"  if $D;
        return [];
    }

    warn "<<<< FAILED >>>>\n"  if $D;
    #no match
    return [];
}

#return a set of HSPs suitable for tiling, that are consistent with the query
#and sbjct sequence numberings.
sub get_ranked_hsps {
    my ($self, $match, $coll, $key, $qorient) = @_;

    my $tmp = $self->get_hsp_groups($match);

    return []  unless keys %$tmp;

    my $alist = $self->get_ranked_hsps_by_query($tmp, $coll, $key, $qorient);

    return $self->combine_frags_by_centroid($alist, $qorient);
}

# #selects the first HSP: minimal, first is assumed to be the best one
# sub get_first_ranked_hsp {
#     my ($self, $match, $coll, $key, $qorient) = @_;
#     my $tmp = $self->get_hsp_groups($match);
#     my $alist = $self->get_ranked_hsps_by_query($tmp, $coll, $key, $qorient);
#     return [ $alist->[0] ]  if @$alist;
#     return $alist;
# }

# #selects all HSPs: will include out of sequence order HSPs
# sub get_all_ranked_hsps {
#     my ($self, $match, $coll, $key, $qorient) = @_;
#     my $tmp = $self->get_hsp_groups($match);
#     my $alist = $self->get_ranked_hsps_by_query($tmp, $coll, $key, $qorient);
#     return $alist;
# }


###########################################################################
###########################################################################
package Bio::MView::Build::Row::BLAST;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub posn1 {
    my $qfm = $_[0]->{'seq'}->fromlabel1;
    my $qto = $_[0]->{'seq'}->tolabel1;
    return "$qfm:$qto";
}

sub posn2 {
    my $hfm = $_[0]->{'seq'}->fromlabel2;
    my $hto = $_[0]->{'seq'}->tolabel2;
    return "$hfm:$hto"  if defined $_[0]->num and $_[0]->num;
    return '';
}

#Sort fragments; called by Row::assemble
#sub sort { $_[0]->sort_none }
#sub sort { $_[0]->sort_worst_to_best }
sub sort { $_[0]->sort_best_to_worst }

# #Don't sort fragments; blast default output lists HSPs by decreasing score,
# #so this would be good, but the mview fragment collecting algorithm called by
# #get_ranked_hsps() changes the order.
# sub sort_none {$_[0]}

# #Sort fragments by (1) increasing score, (2) increasing length; used up to
# #MView version 1.58.1, but inconsistent with NO OVERWRITE policy in
# #Sequence.pm
# sub sort_worst_to_best {
#     $_[0]->{'frag'} = [
#         sort {
#             #warn "$a->[7] <=> $b->[7]\n";
#             my $c = $a->[7] <=> $b->[7];                       #max score
#             return $c  if $c != 0;
#             return length(${$a->[0]}) <=> length(${$b->[0]});  #max length
#         } @{$_[0]->{'frag'}}
#        ];
#     $_[0];
# }

#Sort fragments by (1) decreasing score, (2) decreasing length, (3) leftmost
#start position (arbitrary deal breaker). This is consistent with a best-first
#painting algorithm with NO OVERWRITE in Sequence.pm. Note: query fragment
#score is always set to 1 in the parse_* routines.
sub sort_best_to_worst {

    # warn "sort:\n";
    # foreach my $f (@{$_[0]->{'frag'}}) { warn "@$f\n" }

    $_[0]->{'frag'} = [
        sort {
            my $c;

            #warn "$a->[7] <=> $b->[7] (score)\n";
            $c = $b->[7] <=> $a->[7];              #max score
            return $c  if $c != 0;

            my ($alen, $blen) = (length(${$a->[0]}), length(${$b->[0]}));
            #warn "$alen <=> $blen (length)\n";
            $c = $blen <=> $alen;                  #max length
            return $c  if $c != 0;

            #assume all frags have same orientation - they should have
            if ($a->[1] < $a->[2]) {               #a is forward
                #warn "$a->[1] <=> $b->[1] (start, forward)\n";
                $c = $a->[1] <=> $b->[1];          #min start position
                return $c;
            } else {                               #reversed
                #warn "$b->[1] <=> $a->[1] (start, was reversed)\n";
                $c = $b->[1] <=> $a->[1];          #max start position
                return $c;
            }

        } @{$_[0]->{'frag'}}
        ];

    # warn "sort: (finished)\n";
    # foreach my $f (@{$_[0]->{'frag'}}) { warn "@$f\n" }

    $_[0];
}

#wrapper for logical symmetry with Bio::MView::Build::Row::BLASTX
sub assemble {
    my $self = shift;
    $self->SUPER::assemble(@_);
}


###########################################################################
package Bio::MView::Build::Row::BLASTX;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::BLAST);

#recompute range for translated sequence
sub range {
    my $self = shift;
    my ($lo, $hi) = $self->SUPER::range;
    $self->translate_range($lo, $hi);
}

#assemble translated
sub assemble {
    my $self = shift;
    foreach my $frag (@{$self->{'frag'}}) {
        ($frag->[1], $frag->[2]) =
            $self->translate_range($frag->[1], $frag->[2]);
    }
    $self->SUPER::assemble(@_);
}


###########################################################################
1;
