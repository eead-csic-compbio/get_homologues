# Copyright (C) 1996-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# FASTA 1, 2, 3
#
# The fasta programs can produce different gapped alignments for the same
# database hit. These have the same identifier in the RANK section. Here,
# they are resolved by building a (unique?) key from identifier, initn, and
# init1. This Build doesn't attempt to merge such fragments as done for
# gapped BLAST, ie., each FASTA alignment is taken as a complete solution.
#
###########################################################################
package Bio::MView::Build::Format::FASTA;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Search;
use Bio::Util::Regexp;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Search);

#the name of the underlying Bio::Parse::Format parser
sub parser { 'FASTA' }

my %Known_Parameters =
    (
     #name        => [ format  default ]
     'minopt'     => [ $RX_Ureal,  undef   ],

     #GCG FASTA (version 2)
     'strand'     => [ [],         undef   ],
    );

#our own constructor since this is the entry point for different subtypes
sub new {
    shift;  #discard our own type
    my $self = new Bio::MView::Build::Search(@_);
    my ($type, $p, $v, $file);

    #determine the real type from the underlying parser
    ($p, $v) = (lc $self->{'entry'}->{'format'}, $self->{'entry'}->{'version'});

    #if this is a pre-3.3 fasta call the FASTA2 parser
    if ($v eq '3') {
        my $header = $self->{'entry'}->parse(qw(HEADER));
        $v = 2  if $header->{'version'} =~ /^3\.(\d+)/ and $1 < 3;
    }

    $type = "Bio::MView::Build::Format::FASTA$v";
    ($file = $type) =~ s/::/\//g;
    require "$file.pm";

    $type .= "::$p";
    bless $self, $type;

    $self->initialise;

    $self;
}

#called by the constructor
sub initialise {
    my $self = shift;
    #warn "initialise\n";
    #schedule by strand orientation
    $self->{scheduler} = new Bio::MView::Build::Scheduler([qw(+ -)]);
    $self;
}

#called on each iteration
sub reset_child {
    my $self = shift;
    #warn "reset_child [@{[$PAR->get('strand')]}]\n";
    $self->{scheduler}->filter($PAR->get('strand'));
    $self;
}

#current strand being processed
sub strand { $_[0]->{scheduler}->item }

#compare argument with current strand
sub use_strand { $_[0]->{scheduler}->use_item($_[1]) }

#minopt filter
sub skip_frag {
    my ($self, $opt) = @_;
    return 1  if defined $PAR->get('minopt') and $opt < $PAR->get('minopt');
    return 0;
}

#remove query and hit columns at gaps and frameshifts in the query sequence;
#downcase the bounding hit symbols in the hit sequence thus affected and,
#for frameshifts, downcase the bounding symbols in the query too. remove
#leading/trailing space from the query.
sub strip_query_gaps {
    my ($self, $query, $sbjct, $leader, $trailer) = @_;

    my $gapper = sub {
        my ($query, $sbjct, $char, $doquery) = @_;

        while ( (my $i = index($$query, $char)) >= 0 ) {

            #downcase preceding symbol
            if (defined substr($$query, $i-1, 1)) {
                substr($$query, $i-1, 1) = lc substr($$query, $i-1, 1)
                    if $doquery;
                substr($$sbjct, $i-1, 1) = lc substr($$sbjct, $i-1, 1);
            }

            #consume more of same in query and hit
            while (substr($$query, $i, 1) eq $char) {
                substr($$query, $i, 1) = '';
                substr($$sbjct, $i, 1) = '';
            }

            #downcase succeeding symbol
            if (defined substr($$query, $i, 1)) {
                substr($$query, $i, 1) = lc substr($$query, $i, 1)
                    if $doquery;
                substr($$sbjct, $i, 1) = lc substr($$sbjct, $i, 1);
            }
        }
    };

    #warn "sqg(in  q)=[$$query]\n";
    #warn "sqg(in  h)=[$$sbjct]\n";

    &$gapper($query, $sbjct, '-',  0);  #mark gaps in sbjct only
    &$gapper($query, $sbjct, '/',  1);  #mark frameshifts in both
    &$gapper($query, $sbjct, '\\', 1);  #mark frameshifts in both

    #strip query terminal white space
    $trailer = length($$query) - $leader - $trailer;
    $$query  = substr($$query, $leader, $trailer);
    $$sbjct  = substr($$sbjct, $leader, $trailer);

    #replace sbjct leading/trailing white space with gaps
    $$sbjct =~ s/\s/-/g;

    #warn "sqg(out q)=[$$query]\n";
    #warn "sqg(out h)=[$$sbjct]\n";

    $self;
}

#return a hash of frag sets; each set is a list of ALN with a given
#query/sbjct ordering, score, and significance; these sets have multiple
#alternative keys: "qorient/sorient/opt/sig", "qorient/sorient/init1/sig",
#"qorient/sorient/initn/sig", "qorient/sorient/score/sig" to handle variations
#in the input data.
sub get_frag_groups {
    my ($self, $match) = @_;
    my $hash = {};
    my $sum = $match->parse(qw(SUM));
    foreach my $aln ($match->parse(qw(ALN))) {
        my $sig     = exists $sum->{'expect'} ? $sum->{'expect'} : '-1';
        my $qorient = $aln->{'query_orient'};
        my $sorient = $aln->{'sbjct_orient'};
        my (%tmp, $key);

        if (exists $sum->{'opt'}) {
            $tmp { join("/", $qorient, $sorient, $sum->{'opt'},   $sig) }++;
        }

        if (exists $sum->{'init1'}) {
            $tmp { join("/", $qorient, $sorient, $sum->{'init1'}, $sig) }++;
        }

        if (exists $sum->{'initn'}) {
            $tmp { join("/", $qorient, $sorient, $sum->{'initn'}, $sig) }++;
        }

        if (exists $sum->{'score'}) {
            $tmp { join("/", $qorient, $sorient, $sum->{'score'}, $sig) }++;
        }

        foreach my $key (keys %tmp) {
            push @{ $hash->{$key} }, $aln;
        }
    }
    #warn "SAVE: @{[sort keys %$hash]}\n";
    return $hash;
}

#lookup a frag set in a frag dictionary, trying various keys based on scores
#and significance.
sub get_ranked_frags_by_query {
    my ($self, $hash, $coll, $rkey, $qorient) = @_;

    my $sig = -1;
    if ($coll->item($rkey)->has('expect')) {
        $sig = $coll->item($rkey)->get_val('expect');
    }
    my $qorient2 = ($qorient eq '+' ? '-' : '+');

    my $D = 0;

    my $lookup_o_score_sig = sub {
        my ($qorient, $score, $sig) = @_;
        foreach my $sorient (qw(+ -)) {
            my $key = join("/", $qorient, $sorient, $score, $sig);
            return $hash->{$key}  if exists $hash->{$key};
        }
        return undef;
    };

    my $match;

    warn "KEYS($rkey): @{[sort keys %$hash]}\n"  if $D;

    my $item = $coll->item($rkey);

    foreach my $key (qw(opt init1 initn score)) {
        if ($item->has($key)) {
            my $score = $item->get_val($key);
            #match (score, sig) in query orientation?
            warn "TRY: $qorient $score $sig\n"  if $D;
            $match = &$lookup_o_score_sig($qorient, $score, $sig);
            warn "MATCH (@{[scalar @$match]})\n"  if defined $match and $D;
            return $match  if defined $match;
        }
    }

    warn "<<<< FAILED >>>>\n"  if $D;
    #no match
    return [];
}

#return a set of frags suitable for tiling, that are consistent with the query
#and sbjct sequence numberings.
sub get_ranked_frags {
    my ($self, $match, $coll, $key, $qorient) = @_;
    my $tmp = $self->get_frag_groups($match);

    return []  unless keys %$tmp;

    my $alist = $self->get_ranked_frags_by_query($tmp, $coll, $key, $qorient);

    return $self->combine_frags_by_centroid($alist, $qorient);
}

sub get_scores {
    my ($self, $list, $sum) = @_;

    return unless @$list;

    my ($qorient, $sorient) = ('?', '?');

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
    my $sig = exists $sum->{'expect'} ? $sum->{'expect'} : '';

    ##tweak the significance
    #if ($sig !~ /e/i) {
    #    $sig = sprintf("%.5f", $sig);
    #    $sig =~ s/\.(\d{2}\d*?)0+$/.$1/;  #drop trailing zeros after 2dp
    #    $sig = "0.0"  if $sig == 0;
    #    $sig = "1.0"  if $sig == 1;
    #}

    ($n, $sig, $qorient, $sorient);
}


###########################################################################
###########################################################################
package Bio::MView::Build::Row::FASTA;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub posn1 {
    my $qfm = $_[0]->{'seq'}->fromlabel1;
    my $qto = $_[0]->{'seq'}->tolabel1;
    return "$qfm:$qto"  if defined $qfm and defined $qto;
    return '';
}

sub posn2 {
    my $hfm = $_[0]->{'seq'}->fromlabel2;
    my $hto = $_[0]->{'seq'}->tolabel2;
    return "$hfm:$hto"  if defined $hfm and defined $hto;
    return '';
}

sub assemble {
    my $self = shift;
    $self->SUPER::assemble(@_);
}


###########################################################################
package Bio::MView::Build::Row::FASTX;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA);

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
