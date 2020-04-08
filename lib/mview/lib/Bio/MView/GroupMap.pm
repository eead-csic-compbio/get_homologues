# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::GroupMap;

use Bio::MView::Color::ColorMap;
use Bio::MView::Sequence::Base;
use Exporter;

use vars qw(@ISA @EXPORT $GROUPMAP);

@ISA = qw(Exporter);

@EXPORT = qw($GROUPMAP);

my $GroupMap          = {};   #static hash of consensus schemes
my $Group_Any         = '.';  #key for non-consensus group
my $Default_Group_Any = '.';  #default symbol for non-consensus group

my $Default_PRO_Group = 'P1';
my $Default_DNA_Group = 'D1';

$GROUPMAP = new Bio::MView::GroupMap;  #unique global instance

sub new {
    my $type = shift;
    if (defined $GROUPMAP) {
        die "Bio::MView::GroupMap: instance already exists\n";
    }
    my $self = {};
    bless $self, $type;

    $self->{'map'} = $GroupMap;

    return $GROUPMAP = $self;
}

######################################################################
# public class methods
######################################################################
sub list_groupmap_names { return join(",", sort keys %$GroupMap) }

sub check_groupmap {
    my $map = uc shift;
    return $map  if exists $GroupMap->{$map};
    return undef;
}

sub get_default_groupmap {
    return $Default_PRO_Group  if !defined $_[0];
    return $Default_PRO_Group  if $_[0] eq 'aa';
    return $Default_DNA_Group;
}

sub load_groupmaps {
    my ($stream, $override) = (@_, 1);
    my ($state, $map, $de) = ('start');
    local $_;
    while (<$stream>) {
        my $mapignore;

        #comments, blank lines
        if (/^\s*\#/ or /^\s*$/) {
            next  if $state ne 'map';
            $de .= $_;
            next;
        }

        chomp;

        #group [name]
        if (/^\s*\[\s*(\S+)\s*\]/) {
            $state = 'map';
            $map = uc $1;
            if (exists $GroupMap->{$map} and !$override) {
                $mapignore = 1;  #just for duration of this map
            } else {
                $mapignore = 0;
            }
            $de = '';
            next;
        }

        die "load_groupmaps: groupname undefined\n"  unless defined $map;

        next  if $mapignore;  #forget it if we're not allowing overrides

        #save map description?
        $GroupMap->{$map}->[3] = $de  if $state eq 'map';

        #Group_Any symbol (literal in regexp)
        if (/^\s*\.\s*=>\s*(\S+|\'[^\']+\')/) {
            $state = 'symbol';
            make_group($map, $Group_Any, $1, '') or
                die "load_groupmaps: bad format in line '$_'\n";
            next;
        }

        #general class membership
        if (/^\s*(\S+)\s*=>\s*(\S+|\'[^\']+\')\s*\{\s*(.*)\s*\}/) {
            $state = 'symbol';
            make_group($map, $1, $2, $3) or
                die "load_groupmaps: bad format in line '$_'\n";
            next;
        }

        #trivial class self-membership: different symbol
        if (/^\s*(\S+)\s*=>\s*(\S+|\'[^\']+\')/) {
            $state = 'symbol';
            make_group($map, $1, $2, $1) or
                die "load_groupmaps: bad format in line '$_'\n";
            next;
        }

        #trivial class self-membership: same symbol
        if (/^\s*(\S+)/) {
            $state = 'symbol';
            make_group($map, $1, $1, $1) or
                die "load_groupmaps: bad format in line '$_'\n";
            next;
        }

        #default
        die "load_groupmaps: bad format in line '$_'\n";
    }
    close $stream;

    #set default symbol for each groupmap
    foreach my $map (keys %$GroupMap) {
        make_group($map, $Group_Any, $Default_Group_Any, '')
            unless exists $GroupMap->{$map}->[0]->{$Group_Any};
    }
}

#return a descriptive listing of supplied groups or all groups
sub dump_groupmaps {
    my $html = shift;

    my ($s, $c0, $c1, $c2) = ('', '', '', '');

    ($c0, $c1, $c2) = (
        "<SPAN style=\"color:$Bio::MView::Color::ColorMap::Colour_Black\">",
        "<SPAN style=\"color:$Bio::MView::Color::ColorMap::Colour_Comment\">",
        "</SPAN>")  if $html;

    $s .= "$c1#Consensus group listing - suitable for reloading.\n";
    $s .= "#Character matching is case-insensitive.\n";
    $s .= "#Non-consensus positions default to '$Default_Group_Any' symbol.\n";
    $s .= "#Sequence gaps are shown as ' ' (space) symbols.$c2\n\n";

    @_ = keys %$GroupMap  unless @_;

    foreach my $group (sort @_) {
        $s .= "$c0\[$group]$c2\n";
        $s .= "$c1$GroupMap->{$group}->[3]";
        $s .= "#description =>  symbol  members$c2\n";

        my $p = $GroupMap->{$group}->[0];

        foreach my $class (sort keys %{$p}) {

            next  if $class eq '';  #gap character

            my $sym = $p->{$class}->[0];

            #wildcard
            if ($class eq $Group_Any) {
                $sym = "'$sym'"  if $sym =~ /\s/;
                $s .= sprintf "%-12s =>  %-6s\n", $class, $sym;
                next;
            }

            #consensus symbol
            $sym = "'$sym'"  if $sym =~ /\s/;
            $s .= sprintf "%-12s =>  %-6s  { ", $class, $sym;
            $s .= join(", ", sort keys %{$p->{$class}->[1]}) . " }\n";
        }
        $s .= "\n";
    }

    return $s;
}

#compute consensus group membership for one alignment column; returns
#a hash of consensus group percentages
sub tally_column {
    my ($gname, $col, $gaps) = (@_, 1);

    unless ($GROUPMAP->has_groupmap($gname)) {
        die "tally_column: unknown consensus group '$gname'\n";
    }

    #warn "tally_column: $gname\n";

    my $group = $GROUPMAP->{'map'}->{$gname}->[0];
    my $score = {};
    my $depth;

    #initialise tallies
    foreach my $class (keys %$group) { $score->{$class} = 0 }

    #select score normalization
    if ($gaps) {
        #by total number of rows (sequence + non-sequence)
        $depth = @$col;
    } else {
        #by rows containing sequence in this column
        $depth = 0;
        map { $depth++  if Bio::MView::Sequence::is_char(0, $_) } @$col;
    }
    #warn "($group, [@$col], $gaps, $depth)\n";

    #empty column? use gap symbol
    if ($depth < 1) {
        $score->{''} = 100;
        return $score;
    }

    #tally class scores by column symbol (except gaps), which is upcased
    foreach my $class (keys %$group) {
        foreach my $sym (@$col) {
            next  unless Bio::MView::Sequence::is_char(0, $sym) or $gaps;
            $score->{$class}++  if exists $group->{$class}->[1]->{uc $sym};
        }
        $score->{$class} = 100.0 * $score->{$class} / $depth;
    }

    return $score;
}

#given a columnwise consensus group tally hash, return the consensus
#string for this particular row
sub consensus {
    my ($tally, $gname, $threshold, $ignore) = @_;

    unless ($GROUPMAP->has_groupmap($gname)) {
        die "consensus: unknown consensus group '$gname'\n";
    }

    my $group = $GROUPMAP->{'map'}->{$gname}->[0];
    my $consensus = '';

    #iterate over all columns
    for (my $i=0; $i<@$tally; $i++) {

        my ($score, $bstscore, $bstclass) = ($tally->[$i], 0, undef);

        #iterate over all allowed subsets
        foreach my $class (keys %$group) {

            next  if $class eq $Group_Any;  #wildcard

            if ($class ne '') {
                #non-gap classes: may want to ignore certain classes
                next  if $ignore eq 'singleton' and $class eq $group->{$class}->[0];

                next  if $ignore eq 'class'     and $class ne $group->{$class}->[0];
            }

            #choose smallest class exceeding threshold and highest
            #score percentage when same size

            #warn "[$i] $class, $score->{$class}\n";

            if ($score->{$class} >= $threshold) {

                #first pass
                if (! defined $bstclass) {
                    $bstclass = $class;
                    $bstscore = $score->{$class};
                    next;
                }

                #larger? this set should be rejected
                if (keys %{$group->{$class}->[1]} >
                    keys %{$group->{$bstclass}->[1]}) {
                    next;
                }

                #smaller? this set should be kept
                if (keys %{$group->{$class}->[1]} <
                    keys %{$group->{$bstclass}->[1]}) {
                    $bstclass = $class;
                    $bstscore = $score->{$class};
                    next;
                }

                #same size: new set has better score?
                if ($score->{$class} > $bstscore) {
                    $bstclass = $class;
                    $bstscore = $score->{$class};
                    next;
                }
            }
        }

        if (defined $bstclass) {
            if ($bstclass eq '' and $bstscore < 100) {
                $bstclass = $Group_Any #some non-gaps
            }
        } else {
            $bstclass = $Group_Any;  #wildcard
        }
        #warn "DECIDE [$i] '$bstclass' $bstscore [$group->{$bstclass}->[0]]\n";
        $consensus .= $group->{$bstclass}->[0];
    }

    return \$consensus;
}

######################################################################
# public methods
######################################################################
sub has_groupmap {
    my ($self, $map) = @_;
    return 1  if exists $self->{'map'}->{$map};
    return 0;
}

sub is_singleton {
    my ($self, $map, $cg) = @_;
    return 0  unless exists $self->{'map'}->{$map};
    return 0  unless exists $self->{'map'}->{$map}->[2]->{$cg};
    return 0  unless keys %{ $self->{'map'}->{$map}->[2]->{$cg} } == 1;
    return 1;
}

sub in_consensus_group {
    my ($self, $map, $cg, $c) = @_;
    return 0  unless exists $self->{'map'}->{$map};
    return 0  unless exists $self->{'map'}->{$map}->[1]->{$c};
    return 0  unless exists $self->{'map'}->{$map}->[1]->{$c}->{$cg};
    return 1;
}

######################################################################
# private methods
######################################################################
sub make_group {
    my ($map, $class, $sym, $members) = @_;

    $sym =~ s/^\'//;
    $sym =~ s/\'$//;
    return  0  if length $sym > 1;  #fail

    $members = uc $members;
    $members =~ s/[\s,]//g;
    $members =~ s/''/ /g;
    $members = [ split(//, $members) ];

    #class => symbol
    $GroupMap->{$map}->[0]->{$class}->[0] = $sym;

    foreach my $m (@$members) {
        next  unless defined $m;

        #class  => member existence
        $GroupMap->{$map}->[0]->{$class}->[1]->{$m} = 1;

        #member => symbol existence
        $GroupMap->{$map}->[1]->{$m}->{$sym} = 1;

        #symbol => members
        $GroupMap->{$map}->[2]->{$sym}->{$m} = 1;
    }

    return 1;  #ok
}

######################################################################
# debug
######################################################################
sub dump_group {
    push @_, keys %$GroupMap    unless @_;
    my ($group, $class, $mem, $p);
    warn "Groups by class\n";
    foreach $group (@_) {
        warn "[$group]\n";
        $p = $GroupMap->{$group}->[0];
        foreach $class (keys %{$p}) {
            warn "$class  =>  $p->{$class}->[0]  { ",
                join(" ", keys %{$p->{$class}->[1]}), " }\n";
        }
    }
    warn "Groups by membership\n";
    foreach $group (@_) {
        warn "[$group]\n";
        $p = $GroupMap->{$group}->[1];
        foreach $mem (keys %{$p}) {
            warn "$mem  =>  { ", join(" ", keys %{$p->{$mem}}), " }\n";
        }
    }
}

######################################################################
eval { load_groupmaps(\*DATA) }; warn($@), exit 1  if $@;

1;

__DATA__
#label       =>  symbol  { member list }

[P1]
#Protein consensus: conserved physicochemical classes, derived from
#the Venn diagrams of: Taylor W. R. (1986). The classification of amino acid
#conservation. J. Theor. Biol. 119:205-218.
.            =>  .
G            =>  G       { G }
A            =>  A       { A }
I            =>  I       { I }
V            =>  V       { V }
L            =>  L       { L }
M            =>  M       { M }
F            =>  F       { F }
Y            =>  Y       { Y }
W            =>  W       { W }
H            =>  H       { H }
C            =>  C       { C }
P            =>  P       { P }
K            =>  K       { K }
R            =>  R       { R }
D            =>  D       { D }
E            =>  E       { E }
Q            =>  Q       { Q }
N            =>  N       { N }
S            =>  S       { S }
T            =>  T       { T }
aromatic     =>  a       { F, Y, W, H }
aliphatic    =>  l       { I, V, L }
hydrophobic  =>  h       { I, V, L,   F, Y, W, H,   A, G, M, C, K, R, T }
positive     =>  +       { H, K, R }
negative     =>  -       { D, E }
charged      =>  c       { H, K, R,   D, E }
polar        =>  p       { H, K, R,   D, E,   Q, N, S, T, C }
alcohol      =>  o       { S, T }
tiny         =>  u       { G, A, S }
small        =>  s       { G, A, S,   V, T, D, N, P, C }
turnlike     =>  t       { G, A, S,   H, K, R, D, E, Q, N, T, C }
stop         =>  *       { * }

[D1]
#DNA consensus: conserved ring types
#Ambiguous base R is purine: A or G
#Ambiguous base Y is pyrimidine: C or T or U
.            =>  .
A            =>  A       { A }
C            =>  C       { C }
G            =>  G       { G }
T            =>  T       { T }
U            =>  U       { U }
purine       =>  r       { A, G,  R }
pyrimidine   =>  y       { C, T, U,  Y }

[CYS]
#Protein consensus: conserved cysteines
.            =>  .
C            =>  C       { C }

###########################################################################
