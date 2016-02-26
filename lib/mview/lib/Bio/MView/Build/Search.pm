# Copyright (C) 1997-2006 Nigel P. Brown
# $Id: Search.pm,v 1.12 2015/06/18 21:26:11 npb Exp $

###########################################################################
package Bio::MView::Build::Search;

use Bio::MView::Build;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build);

#there's a query sequence to account for an extra row
sub has_query {1}

sub use_row {
    my ($self, $num, $nid, $sid) = @_;

    #warn "use_row($num, $nid, $sid)\n";
    
    #first, check explicit keeplist and reference row
    foreach my $pat (@{$self->{'keeplist'}}, $self->{'ref_id'}) {

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
    foreach my $pat (@{$self->{'skiplist'}}, $self->{'ref_id'}) {

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

sub map_id {
    my ($self, $ref) = @_;
    $ref = 0  if $ref =~ /query/i;
    $self->SUPER::map_id($ref);
}

#selects the first frag as a seed then add remaining frags that satisfy query
#and sbjct sequence ordering constraints
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

sub key {
    my $self = shift;
    join($DELIM, @_);
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
    $self;
}

sub has {
    my ($self, $key) = @_;
    return 1  if exists $self->{keys}->{$key};
    return 0;
}

sub item {
    my ($self, $key) = @_;
    die "get: unknown key '$key'\n"  unless exists $self->{keys}->{$key};
    $self->{list}->[ $self->{keys}->{$key} ];
}

sub add_frags {
    my ($self, $key, $qf, $qt, $qdata, $hdata) = @_;
    my ($q, $q1, $q2, @qrest) = @$qdata;
    #warn "[$q1, $q2, @qrest]\n";
    $self->item($QUERY)->add_frag($q, $qf,$qt, $q1,$q2, 0,0, @qrest);
    my ($h, $h1, $h2, @hrest) = @$hdata;
    #warn "[$h1, $h2, @hrest]\n";
    $self->item($key)->add_frag($h, $qf,$qt, $q1,$q2, $h1,$h2, @hrest);
    $self;
}

sub list {
    my $self = shift;
    $self->_discard_empty_ranges;
    $self->{list};
}

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
