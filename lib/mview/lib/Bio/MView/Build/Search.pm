# Copyright (C) 1997-2006 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Build::Search;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Base;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Base);

######################################################################
# public methods
######################################################################
#override: the query sequence accounts for an extra row
sub is_search {1}

#selects the first frag as a seed then add remaining frags that satisfy
#query and sbjct sequence ordering constraints
sub combine_frags_by_centroid {
    my ($self, $alist, $qorient) = @_;

    return $alist  if @$alist < 2;  #0 or 1

    my $strictly_ordered = sub {
        my ($o, $a, $b) = @_;
        return $a < $b  if $o eq '+';
        return $a > $b  if $o eq '-';
        return undef;
    };

    my $coincident = sub { $_[0] == $_[1] and $_[2] == $_[3] };

    my $lo = sub {
        my ($o, $a, $b) = @_;
        return ($a < $b ? $a : $b)  if $o eq '+';
        return ($a > $b ? $a : $b)  if $o eq '-';
        return undef;
    };

    my $hi = sub {
        my ($o, $a, $b) = @_;
        return ($a > $b ? $a : $b)  if $o eq '+';
        return ($a < $b ? $a : $b)  if $o eq '-';
        return undef;
    };

    my $query_centre = sub {
        ($_[0]->{'query_start'} + $_[0]->{'query_stop'}) / 2.0;
    };

    my $sbjct_centre = sub {
        ($_[0]->{'sbjct_start'} + $_[0]->{'sbjct_stop'}) / 2.0;
    };

    my $sort_downstream = sub {
        my ($o, $alist) = @_;

        #sort on centre position (increasing) and length (decreasing)
        my $sort = sub {
            my $av = &$query_centre($a);
            my $bv = &$query_centre($b);
            return $av <=> $bv  if $av != $bv;  #centre
            $av = abs($a->{'query_start'} - $a->{'query_stop'});
            $bv = abs($b->{'query_start'} - $b->{'query_stop'});
            return $bv <=> $av;  #length: choose larger, so flip order
        };

        my @tmp = sort $sort @$alist;
        return @tmp          if $o eq '+';  #increasing
        return reverse @tmp  if $o eq '-';  #decreasing
        return @$alist;
    };

    my $sort_upstream = sub { reverse &$sort_downstream(@_) };

    my $aln = shift @$alist;  #first element with the best score
    my ($qm, $sm) = (&$query_centre($aln), &$sbjct_centre($aln));
    my $sorient = $aln->{'sbjct_orient'};
    my @tmp = ($aln);

    #my $D = 0;
    #warn "\n"  if $D;

    #accept the first frag then process the remainder rejecting any whose
    #query or sbjct midpoint does not extend the alignment monotonically;
    #effectively, we grow the alignment diagonally out from boxes at each end
    #of the growing chain in the 2D alignment matrix.

    #grow upstream
    foreach my $aln (&$sort_upstream($qorient, $alist)) {

        my ($aqm, $asm) = (&$query_centre($aln), &$sbjct_centre($aln));

        #warn "UP($qorient$sorient): centres $aqm .. $qm : $asm .. $sm  @{[&$strictly_ordered($qorient, $aqm, $qm)?'1':'0']} @{[&$strictly_ordered($sorient, $asm, $sm)?'1':'0']}\n"  if $D;

        next  unless
            (&$strictly_ordered($qorient, $aqm, $qm) and
             &$strictly_ordered($sorient, $asm, $sm))
            or
            &$coincident($aqm, $qm, $asm, $sm);

        ($qm, $sm) = (&$lo($qorient, $aqm, $qm), &$lo($sorient, $asm, $sm));
        #warn "up($qorient$sorient): [$aln->{'sbjct'}]\n"  if $D;
        push @tmp, $aln;
    }

    #grow downstream
    foreach my $aln (&$sort_downstream($qorient, $alist)) {

        next  if grep {$_ eq $aln} @tmp;  #seen it already

        my ($aqm, $asm) = (&$query_centre($aln), &$sbjct_centre($aln));

        #warn "DN($qorient$sorient): centres $qm .. $aqm : $sm .. $asm  @{[&$strictly_ordered($qorient, $qm, $aqm)?'1':'0']} @{[&$strictly_ordered($sorient, $sm, $asm)?'1':'0']}\n"  if $D;

        next  unless
            (&$strictly_ordered($qorient, $qm, $aqm) and
             &$strictly_ordered($sorient, $sm, $asm))
            or
            &$coincident($aqm, $qm, $asm, $sm);

        ($qm, $sm) = (&$hi($qorient, $aqm, $qm), &$hi($sorient, $asm, $sm));
        #warn "dn($qorient$sorient): [$aln->{'sbjct'}]\n"  if $D;
        push @tmp, $aln;
    }
    return \@tmp;
}

######################################################################
# protected methods
######################################################################
#override: are sequence inserts being kept?
sub test_if_aligned {
    my ($self, $lo, $hi) = @_;
    return 0  if $PAR->get('keepinserts');
    return 1;
}

#override
sub use_row {
    my ($self, $num, $nid, $sid) = @_;

    #warn "use_row($num, $nid, $sid)\n";

    #first, check explicit keeplist and reference row
    foreach my $pat (@{$PAR->get('keeplist')}, $PAR->get('ref_id')) {

        #Search subclass only
        return 1  if $pat eq '0'     and $num == 0;
        return 1  if $pat eq 'query' and $num == 0;

        #look at row number
        return 1  if $nid eq $pat;      #major OR major.minor
        if ($nid =~ /^\d+$/ and $pat =~ /^(\d+)\./) {
            #major matched by major.*
            return 1  if $nid eq $1;
        } elsif ($pat =~ /^\d+$/ and $nid =~ /^(\d+)\./) {
            #major.* matched by major
            return 1  if $1 eq $pat;
        }

        #look at identifier
        return 1  if $sid eq $pat;      #exact match
        if ($pat =~ /^\/(.*)\/$/) {     #regex match (case insensitive)
            return 1  if $sid =~ /$1/i;
        }
    }

    #second, check skiplist and reference row
    foreach my $pat (@{$PAR->get('skiplist')}, $PAR->get('ref_id')) {

        #Search subclass only
        return 0  if $pat eq '0'     and $num == 0;
        return 0  if $pat eq 'query' and $num == 0;

        #look at row number
        return 0  if $nid eq $pat;      #major OR major.minor
        if ($nid =~ /^\d+$/ and $pat =~ /^(\d+)\./) {
            #major matched by major.*
            return 0  if $nid eq $1;
        } elsif ($pat =~ /^\d+$/ and $nid =~ /^(\d+)\./) {
            #major.* matched by major
            return 0  if $1 eq $pat;
        }

        #look at identifier
        return 0  if $sid eq $pat;      #exact match
        if ($pat =~ /^\/(.*)\/$/) {     #regex match (case insensitive)
            return 0  if $sid =~ /$1/i;
        }
    }

    #assume implicit membership of keeplist
    return 1;    #default
}

#override
sub map_id {
    my ($self, $id) = @_;
    $id = 0  if $id =~ /query/i;
    return $self->SUPER::map_id($id);
}

#subclass overrides: remove query and hit columns at gaps in the query
#sequence and downcase the bounding hit symbols in the hit sequence thus
#affected.
sub strip_query_gaps {
    my ($self, $query, $sbjct) = @_;

    #warn "sqg(in  q)=[$$query]\n";
    #warn "sqg(in  h)=[$$sbjct]\n";

    #no gaps in query
    return    if index($$query, '-') < 0;

    #iterate over query frag symbols
    while ( (my $i = index($$query, '-')) >= 0 ) {

        #downcase preceding symbol in hit
        if (defined substr($$query, $i-1, 1)) {
            substr($$sbjct, $i-1, 1) = lc substr($$sbjct, $i-1, 1);
        }

        #consume gap symbols in query and hit
        while (substr($$query, $i, 1) eq '-') {
            substr($$query, $i, 1) = "";
            substr($$sbjct, $i, 1) = "";
        }

        #downcase succeeding symbol in hit
        if (defined substr($$query, $i, 1)) {
            substr($$sbjct, $i, 1) = lc substr($$sbjct, $i, 1);
        }

        #warn "sqg(out q)=[$$query]\n";
        #warn "sqg(out h)=[$$sbjct]\n";
    }
    $self;
}

###########################################################################
package Bio::MView::Build::Search::Collector;

my $QUERY = 0;
my $DELIM = '.';

sub new {
    my $type = shift;
    my $self = {};
    bless $self, $type;

    $self->{parent} = shift;
    $self->{keys} = {};
    $self->{list} = [];

    $self;
}

######################################################################
# public methods
######################################################################
sub has {
    my ($self, $key) = @_;
    return 1  if exists $self->{keys}->{$key};
    return 0;
}

sub key {
    my $self = shift;
    return join($DELIM, @_);
}

sub item {
    my ($self, $key) = @_;
    die "get: unknown key '$key'\n"  unless exists $self->{keys}->{$key};
    return $self->{list}->[ $self->{keys}->{$key} ];
}

sub list {
    my $self = shift;
    $self->_discard_empty_ranges;
    return $self->{list};
}

sub insert {
    my ($self, $item) = (shift, shift);
    push @_, $QUERY  unless @_;
    foreach my $key (@_) {
        #allow this
        #warn "insert: key already exists '$key'\n"
        #    if exists $self->{keys}->{$key};
        $self->{keys}->{$key} = scalar @{$self->{list}};
    }
    push @{$self->{list}}, $item;
}

sub add_frags {
    my ($self, $key, $qf, $qt, $qdata, $hdata) = @_;
    my ($q, $q1, $q2, @qrest) = @$qdata;
    #warn "[$q1, $q2, @qrest]\n";
    $self->item($QUERY)->add_frag($q, $qf,$qt, $q1,$q2, 0,0, @qrest);
    my ($h, $h1, $h2, @hrest) = @$hdata;
    #warn "[$h1, $h2, @hrest]\n";
    $self->item($key)->add_frag($h, $qf,$qt, $q1,$q2, $h1,$h2, @hrest);
}

######################################################################
# private methods
######################################################################
#remove hits lacking positional data; finally remove the query itself if
#that's all that's left
sub _discard_empty_ranges {
    my $self = shift;
    my $hit = $self->{list};
    for (my $i=1; $i<@$hit; $i++) {
        #warn "hit[$i]= $hit->[$i]->{'cid'} [", scalar @{$hit->[$i]->{'frag'}},"]\n";
        if (@{$hit->[$i]->{'frag'}} < 1) {
            splice(@$hit, $i--, 1);
        }
    }
    pop @$hit  unless @$hit > 1;
    $self;
}

######################################################################
# debug
######################################################################
sub dump {
    my $self = shift;
    warn "Collector:\n";
    warn "keys: @{[scalar keys %{$self->{keys}}]} list: @{[scalar @{$self->{list}}]}\n";
    foreach my $k (sort { $self->{keys}->{$a} <=> $self->{keys}->{$b} } keys %{$self->{keys}}) {
        my $i = $self->{keys}->{$k};
        warn "[$self->{keys}->{$k}] <= $k => $self->{list}->[$i]\n";
    }
}

###########################################################################
1;
