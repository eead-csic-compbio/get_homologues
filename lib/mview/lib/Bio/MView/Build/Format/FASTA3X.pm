# Copyright (C) 1996-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# FASTA 3X (34/35/36)
#
#   fasta, fastx, fasty, tfasta, tfastx, tfasty, tfastxy,
#   fastm, fastf, fasts, ggsearch, glsearch, ssearch
#
###########################################################################
use Bio::MView::Build::Format::FASTA3;

use strict;

###########################################################################
###########################################################################
package Bio::MView::Build::Row::FASTA3X;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::fasta);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::fastx);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::fasty);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfasta);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfastx);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfasty);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfastxy;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3::tfastxy);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fastm;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 1,   1,    'initn',         'initn',      '5N',      ''  ],
    [ 2,   2,    'init1',         'init1',      '5N',      ''  ],
    [ 3,   3,    'bits',          'bits',       '7N',      ''  ],
    [ 4,   4,    'expect',        'E-value',    '9N',      ''  ],
    [ 5,   5,    'sn',            'sn',         '3N',      ''  ],
    [ 6,   6,    'sl',            'sl',         '3N',      ''  ],
    [ 7,   7,    'query_orient',  'qy',         '2S',      '?' ],
    [ 8,   8,    'sbjct_orient',  'ht',         '2S',      '?' ],
    ]
}


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fastf;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::fastm);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::fasts;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::fastm);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfastm;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::fastm);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfastf;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::fastm);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::tfasts;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::fastm);


###########################################################################
package Bio::MView::Build::Row::FASTA3X::ssearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 1,   1,    'score',         'S-W',        '5N',      ''  ],
    [ 2,   2,    'bits',          'bits',       '7N',      ''  ],
    [ 3,   3,    'expect',        'E-value',    '9N',      ''  ],
    [ 4,   4,    'query_orient',  'qy',         '2S',      '?' ],
    [ 5,   5,    'sbjct_orient',  'ht',         '2S',      '?' ],
    ]
}


###########################################################################
package Bio::MView::Build::Row::FASTA3X::ggsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::ssearch);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 1,   1,    'score',         'N-W',        '5N',      ''  ],
    [ 2,   2,    'bits',          'bits',       '7N',      ''  ],
    [ 3,   3,    'expect',        'E-value',    '9N',      ''  ],
    [ 4,   4,    'query_orient',  'qy',         '2S',      '?' ],
    [ 5,   5,    'sbjct_orient',  'ht',         '2S',      '?' ],
    ]
}


###########################################################################
package Bio::MView::Build::Row::FASTA3X::glsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::FASTA3X::ggsearch);


###########################################################################
###########################################################################
package Bio::MView::Build::Format::FASTA3X;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::fasta);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::fastx);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::fasty);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfasta;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfasta);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfastx);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfasty;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfasty);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfastxy;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3::tfastxy);


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fastm;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X); #note

my $FASTMFS_SPACER = '\001' x 5;

#called by the constructor
sub initialise {
    my $self = shift;

    #schedule by query peptide tuple
    my $peplist = $self->parse_query_tuples;

    $self->{scheduler} = new Bio::MView::Build::Scheduler($peplist);
    $self;
}

#called on each iteration
sub reset_child {
    my $self = shift;
    #warn "reset_child [all]\n";
    $self->{scheduler}->filter;
    $self;
}

#current peptide tuple being processed
sub pepnum { $_[0]->{scheduler}->itemnum }
sub peptup { $_[0]->{scheduler}->item }

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Peptide tuple: @{[$self->pepnum]}  @{[$self->peptup]}\n";
    $s;
}

sub makepeptup {
    $_[0] =~ s/^\s+//o;
    $_[0] =~ s/\s+$//o;
    $_[0] =~ s/[-\s]+/, /og;
    $_[0];
}

sub parse_query_tuples {
    my $self = shift;
    my ($peplist, $pephash) = ([], {});

    foreach my $match ($self->{'entry'}->parse(qw(MATCH))) {

        #tfast[mfs] ALN needs to access its SUM
        my $sum = $match->parse(qw(SUM));

        foreach my $aln ($match->parse(qw(ALN))) {

            my $peptup = makepeptup("$aln->{'query'}");

            next  if exists $pephash->{$peptup};

            push @$peplist, $peptup;
            $pephash->{$peptup}++;

            #warn "PEP: $peptup\n";
        }
    }

    $self->{'entry'}->free_parsers(qw(MATCH));

    $peplist;
}

sub parse {
    my $self = shift;

    #warn "\nfasta3X::parse: entering\n";

    #all peptide tuples done?
    if (! defined $self->{scheduler}->next) {
        $self->{'entry'}->free_parsers();
        #warn "fasta3X::parse: exiting at end\n";
        return;
    }

    my $ranking = $self->{'entry'}->parse(qw(RANK));

    #no hits
    if (! defined $ranking) {
        #warn "fasta3X::parse: exiting on empty parse\n";
        return;
    }

    #identify the query
    my $header = $self->{'entry'}->parse(qw(HEADER));

    my $query = 'Query';
    if ($header->{'query'} ne '') {
        $query = $header->{'query'};
    } elsif ($header->{'queryfile'} =~ m,.*/([^\.]+)\.,) {
        $query = $1;
    } else {

    }

    my $coll = new Bio::MView::Build::Search::Collector($self);

    my $rtype = $1  if ref($self) =~ /::([^:]+)$/;
    my $class = "Bio::MView::Build::Row::FASTA3X::$rtype";

    $coll->insert((new $class(
                       '',
                       $query,
                       '',
                       {
                           'initn'        => '',
                           'init1'        => '',
                           'bits'         => '',
                           'expect'       => '',
                           'sn'           => '',
                           'sl'           => '',
                           'query_orient' => '+',
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

        #warn "ADD: [$key1]\n";

        $coll->insert((new $class(
                           $rank,
                           $hit->{'id'},
                           $hit->{'desc'},
                           {
                               'initn'        => $hit->{'initn'},
                               'init1'        => $hit->{'init1'},
                               'bits'         => $hit->{'bits'},
                               'expect'       => $hit->{'expect'},
                               'sn'           => $hit->{'sn'},
                               'sl'           => $hit->{'sl'},
                               'query_orient' => '+',
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

        #ignore hit?
        next  unless $coll->has($key1);

        #warn "USE: [$key1]\n";

        my $aset = $self->get_ranked_frags($match, $coll, $key1, '+');

        #nothing matched
        next  unless @$aset;

        foreach my $aln (@$aset) {
            #apply score/significance filter
            next  if $self->skip_frag($sum->{'initn'});

            #warn "SEE: [$key1]\n";

            #$aln->dump;

            my $peptup = makepeptup("$aln->{'query'}");

            #ignore other peptide tuples
            next  unless $self->{scheduler}->use_item($peptup);

            #for FASTA gapped alignments
            $self->strip_query_gaps(\$aln->{'query'}, \$aln->{'sbjct'},
                                    $aln->{'query_leader'},
                                    $aln->{'query_trailer'});

            my $qstop = $aln->{'query_start'} + length($aln->{'query'}) - 1;

            $coll->add_frags(
                $key1, $aln->{'query_start'}, $qstop, [
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
        #$coll->item($key1)->set_val('expect', $sig); #not used yet
        $coll->item($key1)->set_val('sbjct_orient', $sorient);
    }

    #warn "fasta3X::parse: returning with data\n";

    return $coll->list;
}

#overrides FASTA::strip_query_gaps
sub strip_query_gaps {
    my ($self, $query, $sbjct, $leader, $trailer) = @_;

    my $gapper = sub {
        my ($query, $sbjct, $char) = @_;

        while ( (my $i = index($$query, $char)) >= 0 ) {

            my $pos = $i;

            #downcase preceding symbol
            if (defined substr($$query, $i-1, 1)) {
                substr($$sbjct, $i-1, 1) = lc substr($$sbjct, $i-1, 1);
            }

            #consume more of same in query and hit
            while (substr($$query, $i, 1) eq $char) {
                substr($$query, $i, 1) = '';
                substr($$sbjct, $i, 1) = '';
            }

            #downcase succeeding symbol
            if (defined substr($$query, $i, 1)) {
                substr($$sbjct, $i, 1) = lc substr($$sbjct, $i, 1);
            }

            #insert fixed spacer
            substr($$query, $pos, 0) = $FASTMFS_SPACER;
            substr($$sbjct, $pos, 0) = $FASTMFS_SPACER;
        }
    };

    &$gapper($query, $sbjct, '-');  #mark gaps in sbjct only

    #strip query terminal white space
    $trailer = length($$query) - $leader - $trailer;
    $$query  = substr($$query, $leader, $trailer);
    $$sbjct  = substr($$sbjct, $leader, $trailer);

    #replace sbjct leading/trailing white space with gaps
    $$sbjct =~ s/\s/-/g;

    #replace spacer with '-'
    $$query =~ s/\\001/-/g;
    $$sbjct =~ s/\\001/-/g;

    #warn "sqg(out q)=[$$query]\n";
    #warn "sqg(out h)=[$$sbjct]\n";

    $self;
}


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fastf;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::fastm); #note


###########################################################################
package Bio::MView::Build::Format::FASTA3X::fasts;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::fastm); #note


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfastm;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::fastm); #note


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfastf;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::fastm); #note


###########################################################################
package Bio::MView::Build::Format::FASTA3X::tfasts;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::fastm); #note


###########################################################################
package Bio::MView::Build::Format::FASTA3X::ssearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X); #note

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s  if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Query orientation: " . $self->strand . "\n";
    $s;
}

sub parse {
    my $self = shift;

    #warn "\nfasta3X::parse: entering\n";

    #all strands done?
    if (! defined $self->{scheduler}->next) {
        $self->{'entry'}->free_parsers();
        #warn "fasta3X::parse: exiting at end\n";
        return;
    }

    my $ranking = $self->{'entry'}->parse(qw(RANK));

    #no hits
    if (! defined $ranking) {
        #warn "fasta3X::parse: exiting on empty parse\n";
        return;
    }

    #identify the query
    my $header = $self->{'entry'}->parse(qw(HEADER));

    my $query = 'Query';
    if ($header->{'query'} ne '') {
        $query = $header->{'query'};
    } elsif ($header->{'queryfile'} =~ m,.*/([^\.]+)\.,) {
        $query = $1;
    } else {
        $query = 'Query';
    }

    my $coll = new Bio::MView::Build::Search::Collector($self);

    my $rtype = $1  if ref($self) =~ /::([^:]+)$/;
    my $class = "Bio::MView::Build::Row::FASTA3X::$rtype";

    $coll->insert((new $class(
                       '',
                       $query,
                       '',
                       {
                           'score'        => '',
                           'bits'         => '',
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

        #warn "ADD: [$key1]\n";

        $coll->insert((new $class(
                           $rank,
                           $hit->{'id'},
                           $hit->{'desc'},
                           {
                               'score'        => $hit->{'score'},
                               'bits'         => $hit->{'bits'},
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
            next  if $self->skip_frag($sum->{'score'});

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
        #$coll->item($key1)->set_val('expect', $sig); #not used yet
        $coll->item($key1)->set_val('sbjct_orient', $sorient);
    }

    #warn "fasta3X::parse: returning with data\n";

    return $coll->list;
}


###########################################################################
package Bio::MView::Build::Format::FASTA3X::ggsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::ssearch); #note


###########################################################################
package Bio::MView::Build::Format::FASTA3X::glsearch;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::FASTA3X::ssearch); #note


###########################################################################
1;
