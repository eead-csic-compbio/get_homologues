# Copyright (C) 1997-2013 Nigel P. Brown
# $Id: Consensus.pm,v 1.26 2015/06/14 17:09:03 npb Exp $

###########################################################################
package Bio::MView::Align::Consensus;

use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Align::Row;

use strict;
use vars qw(@ISA
	    $Default_PRO_Colormap $Default_DNA_Colormap $Default_Colormap
	    $Default_PRO_Group $Default_DNA_Group
	    $Default_Group_Any $Default_Ignore $Group);

@ISA = qw(Bio::MView::Align::Sequence);

$Default_PRO_Colormap = 'PC1';    #default colormap name
$Default_DNA_Colormap = 'DC1';    #default colormap name
$Default_Colormap     = $Default_PRO_Colormap;

$Default_PRO_Group    = 'P1';     #default consensus scheme name
$Default_DNA_Group    = 'D1';     #default consensus scheme name
  				  
$Default_Group_Any    = '.';      #default symbol for non-consensus
$Default_Ignore       = '';       #default ignore classes setting
$Group                = {};       #static hash of consensus schemes


my %Known_Ignore_Class = 
    (
     #name
     'none'       => 1,    #don't ignore
     'singleton'  => 1,    #ignore singleton, ie., self-only consensi
     'class'      => 1,    #ignore non-singleton, ie., class consensi
    );

#static load the $Group hash.
eval { load_groupmaps() }; if ($@) {$::COMPILE_ERROR=1; warn $@}


sub load_groupmaps {
    my ($stream, $override) = (@_, \*DATA, 1);
    my ($state, $map, $class, $sym, $members, $c, $de, $mapignore) = (0, {}, undef);
    local $_;
    while (<$stream>) {

	#comments, blank lines
	if (/^\s*\#/ or /^\s*$/) {
	    next  if $state != 1;
	    $de .= $_;
	    next;
	}

	#group [name]
	if (/^\s*\[\s*(\S+)\s*\]/) {
	    $map = uc $1;
	    if (exists $Group->{$map} and !$override) {
		$mapignore = 1;  #just for duration of this map
	    } else {
		$mapignore = 0;
	    }
	    $state = 1;
	    $de = '';
	    next;
	}

	die "Bio::MView::Align::Row::load_groupmaps() groupname undefined\n"    
	    unless defined $map;

	next  if $mapignore;    #forget it if we're not allowing overrides

	#save map description?
	$Group->{$map}->[3] = $de  if $state == 1;

	#ANY symbol
	if (/^\s*\*\s*=>\s*(\S+|\'[^\']+\')/) {
	    $state = 2;
	    $sym     = $1;
	    $sym     =~ s/^\'//;
	    $sym     =~ s/\'$//;
	    chomp; die "Bio::MView::Align::Row::load_groupmaps() bad format in line '$_'\n"    if length $sym > 1;
	    make_group($map, '*', $sym, []);
	}
	
	#general class membership
	if (/^\s*(\S+)\s*=>\s*(\S+|\'[^\']+\')\s*\{\s*(.*)\s*\}/) {
	    $state = 2;
	    ($class, $sym, $members) = ($1, $2, $3);
	    $sym     =~ s/^\'//;
	    $sym     =~ s/\'$//;
	    chomp; die "Bio::MView::Align::Row::load_groupmaps() bad format in line '$_'\n"    if length $sym > 1;
	    $members =~ s/[\s,]//g;
	    $members =~ s/''/ /g;
	    $members = uc $members;
	    $members = [ split(//, $members) ];

	    make_group($map, $class, $sym, $members);

	    next;
	}

	#trivial class self-membership: different symbol
	if (/^\s*(\S+)\s*=>\s*(\S+|\'[^\']+\')/) {
	    $state = 2;
	    ($class, $sym, $members) = ($1, $2, $1);
	    chomp; die "Bio::MView::Align::Row::load_groupmaps() bad format in line '$_'\n"    if length $sym > 1;
	    $members = uc $members;
	    $members = [ split(//, $members) ];

	    make_group($map, $class, $sym, $members);

	    next;
	}

	#trivial class self-membership: same symbol
	if (/^\s*(\S+)/) {
	    $state = 2;
	    ($class, $sym, $members) = ($1, $1, $1);
	    $members = uc $members;
	    $members = [ split(//, $members) ];

	    make_group($map, $class, $sym, $members);

	    next;
	}

	#default
	chomp; die "Bio::MView::Align::Row::load_groupmaps() bad format in line '$_'\n";	
    }
    close $stream;

    foreach $map (keys %$Group) {
	make_group($map, '*', $Default_Group_Any, [])
	    unless exists $Group->{$map}->[0]->{'*'};
	make_group($map, '', $Bio::MView::Sequence::Mark_Spc,
		   [
		    $Bio::MView::Sequence::Mark_Pad, 
		    $Bio::MView::Sequence::Mark_Gap,
		   ]);
    }
}

sub make_group {
    my ($group, $class, $sym, $members) = @_;
    local $_;

    #class => symbol
    $Group->{$group}->[0]->{$class}->[0] = $sym;
    
    foreach (@$members) {
	#class  => member existence
	$Group->{$group}->[0]->{$class}->[1]->{$_} = 1;
	#member => symbol existence
	$Group->{$group}->[1]->{$_}->{$sym} = 1;
	#symbol => members
	$Group->{$group}->[2]->{$sym}->{$_} = 1;
    }
}

sub dump_group {
    my ($group, $class, $mem, $p);
    push @_, keys %$Group    unless @_;
    warn "Groups by class\n";
    foreach $group (@_) {
	warn "[$group]\n";
	$p = $Group->{$group}->[0];
	foreach $class (keys %{$p}) {
	    warn "$class  =>  $p->{$class}->[0]  { ",
		join(" ", keys %{$p->{$class}->[1]}), " }\n";
	}
    }
    warn "Groups by membership\n";
    foreach $group (@_) {
	warn "[$group]\n";
	$p = $Group->{$group}->[1];
	foreach $mem (keys %{$p}) {
	    warn "$mem  =>  { ", join(" ", keys %{$p->{$mem}}), " }\n";
	}
    }
}

#return a descriptive listing of supplied groups or all groups
sub list_groupmaps {
    my $html = shift;
    my ($group, $class, $p, $sym);
    my ($s, $c0, $c1, $c2) = ('', '', '', '');

    ($c0, $c1, $c2) = (
	"<SPAN style=\"color:$Bio::MView::Align::Row::Colour_Black\">",
	"<SPAN style=\"color:$Bio::MView::Align::Row::Colour_Comment\">",
	"</SPAN>")  if $html;

    $s .= "$c1#Consensus group listing - suitable for reloading.\n";
    $s .= "#Character matching is case-insensitive.\n";
    $s .= "#Non-consensus positions default to '$Default_Group_Any' symbol.\n";
    $s .= "#Sequence gaps are shown as ' ' (space) symbols.$c2\n\n";
    
    @_ = keys %$Group  unless @_;

    foreach $group (sort @_) {
	$s .= "$c0\[$group]$c2\n";
	$s .= "$c1$Group->{$group}->[3]";
        $s .= "#description =>  symbol  members$c2\n";

	$p = $Group->{$group}->[0];
	foreach $class (sort keys %{$p}) {
	    
	    next    if $class eq '';    #gap character

	    #wildcard
	    if ($class eq '*') {
		$sym = $p->{$class}->[0];
		$sym = "'$sym'"    if $sym =~ /\s/;
		$s .= sprintf "%-12s =>  %-6s\n", $class, $sym;
		next;
	    }
	    
	    #consensus symbol
	    $sym = $p->{$class}->[0];
	    $sym = "'$sym'"    if $sym =~ /\s/;
	    $s .= sprintf "%-12s =>  %-6s  { ", $class, $sym;
	    $s .= join(", ", sort keys %{$p->{$class}->[1]}) . " }\n";
	}
	$s .= "\n";
    }
    $s;
}

sub check_groupmap {
    if (defined $_[0]) {
	if (exists $Group->{uc $_[0]}) {
	    return uc $_[0];
	}
	return undef;
    }
    return sort keys %$Group;
}

sub check_ignore_class {
    if (defined $_[0]) {
	if (exists $Known_Ignore_Class{lc $_[0]}) {
	    return lc $_[0];
	}
	return undef;
    }
    return map { lc $_ } sort keys %Known_Ignore_Class;
}

sub get_default_groupmap {
    if (! defined $_[0] or $_[0] eq 'aa') {
	#default to protein
	return $Default_PRO_Group;
    }
    #otherwise DNA/RNA explicitly requested
    return $Default_DNA_Group;
}

sub get_color_identity { my $self = shift; $self->SUPER::get_color(@_) }

#colours a row of consensus sequence.
#philosophy:
#  1. give consensus symbols their own colour.
#  2. the consensus may be a residue name: use prevailing residue colour.
#  3. use the prevailing wildcard residue colour.
#  4. give up.
sub get_color_type {
    my ($self, $c, $mapS, $mapG) = @_;
    my ($index, $color, $trans);

    #warn "get_color_type($self, $c, $mapS, $mapG)\n";

    #look in group colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapG}->{$c}) {

	#set transparent(T)/solid(S)
	$trans = $Bio::MView::Align::Colormaps->{$mapG}->{$c}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapG}->{$c}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $mapG\{$c} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");

    }

    #look in sequence colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapS}->{$c}) {

	#set transparent(T)/solid(S)
	$trans = $Bio::MView::Align::Colormaps->{$mapS}->{$c}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapS}->{$c}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $mapS\{$c} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }
    
    #look for wildcard in sequence colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapS}->{'*'}) {

	#set transparent(T)/solid(S)
	$trans = $Bio::MView::Align::Colormaps->{$mapS}->{'*'}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapS}->{'*'}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $mapS\{'*'} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }

    return 0;    #no match
}

#colours a row of 'normal' sequence only where there is a consensus symbol.
#philosophy:
#  1. give residues their own colour.
#  2. use the prevailing wildcard residue colour.
#  3. give up.
sub get_color_consensus_sequence {
    my ($self, $cs, $cg, $mapS, $mapG) = @_;
    my ($index, $color, $trans);

    #warn "get_color_consensus_sequence($self, $cs, $cg, $mapS, $mapG)\n";

    #lookup sequence symbol in sequence colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapS}->{$cs}) {

	#set transparent(T)/solid(S)
	$trans = $Bio::MView::Align::Colormaps->{$mapS}->{$cs}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapS}->{$cs}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$cs/$cg $mapS\{$cs} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }
    
    #lookup wildcard in sequence colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapS}->{'*'}) {

	#set transparent(T)/solid(S)
	$trans = $Bio::MView::Align::Colormaps->{$mapS}->{'*'}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapS}->{'*'}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$cs/$cg $mapS\{'*'} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }

    return 0;    #no match
}

#colours a row of 'normal' sequence using colour of consensus symbol.
#philosophy:
#  1. give residues the colour of the consensus symbol.
#  2. the consensus may be a residue name: use prevailing residue colour.
#  3. use the prevailing wildcard residue colour.
#  4. give up.
sub get_color_consensus_group {
    my ($self, $cs, $cg, $mapS, $mapG) = @_;
    my ($index, $color, $trans);

    #warn "get_color_consensus_group($self, $cs, $cg, $mapS, $mapG)\n";

    #lookup group symbol in group colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapG}->{$cg}) {

	#set transparent(T)/solid(S)/color from GROUP colormap
	$trans = $Bio::MView::Align::Colormaps->{$mapG}->{$cg}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapG}->{$cg}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];
	#warn "$cs/$cg $mapG\{$cg} [$index] [$color] [$trans]\n";

	return ($color, "$trans$index");
    }

    #lookup group symbol in sequence colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapS}->{$cg}) {

	#set transparent(T)/solid(S)/color from SEQUENCE colormap
	$trans = $Bio::MView::Align::Colormaps->{$mapS}->{$cg}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapS}->{$cg}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];
	#warn "$cs/$cg $mapS\{$cg} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }

    #lookup wildcard in SEQUENCE colormap
    if (exists $Bio::MView::Align::Colormaps->{$mapS}->{'*'}) {

	#set transparent(T)/solid(S)/color from SEQUENCE colormap
	$trans = $Bio::MView::Align::Colormaps->{$mapS}->{'*'}->[1];
	$index = $Bio::MView::Align::Colormaps->{$mapS}->{'*'}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];
	#warn "$cs/$cg $mapS\{'*'} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }

    return 0;    #no match
}

sub tally {
    my ($group, $col, $gaps) = (@_, 1);
    my ($score, $class, $sym, $depth) = ({});
    
    if (! exists $Group->{$group}) {
	die "Bio::MView::Align::Consensus::tally() unknown consensus set\n";
    }
    
    #warn "tally: $group\n";

    $group = $Group->{$group}->[0];
    
    #initialise tallies
    foreach $class (keys %$group) { $score->{$class} = 0 }
    
    #select score normalization
    if ($gaps) {
	#by total number of rows (sequence + non-sequence)
	$depth = @$col;
    } else {
	#by rows containing sequence in this column
	$depth = 0;
	map { $depth++ if Bio::MView::Sequence::is_char(0, $_) } @$col;
    }
    #warn "($group, [@$col], $gaps, $depth)\n";

    #empty column? use gap symbol
    if ($depth < 1) {
	$score->{''} = 100;
	return $score;
    }

    #tally class scores by column symbol (except gaps), which is upcased
    foreach $class (keys %$group) {
	foreach $sym (@$col) {
	    next    unless Bio::MView::Sequence::is_char(0, $sym) or $gaps;
	    $score->{$class}++    if exists $group->{$class}->[1]->{uc $sym};
	}
	$score->{$class} = 100.0 * $score->{$class} / $depth;
    }
    $score;
}

sub consensus {
    my ($tally, $group, $threshold, $ignore) = @_;
    my ($class, $bstclass, $bstscore, $consensus, $i, $score);

    if (! exists $Group->{$group}) {
	die "Bio::MView::Align::Consensus::consensus() unknown consensus set\n";
    }
    
    $group = $Group->{$group}->[0];

    $consensus = '';

    #iterate over all columns
    for ($i=0; $i<@$tally; $i++) {
	
	($score, $class, $bstclass, $bstscore) = ($tally->[$i], "", undef, 0);
	
	#iterate over all allowed subsets
	foreach $class (keys %$group) {

	    next    if $class eq '*'; #wildcard
	    
	    if ($class ne '') {
		#non-gap classes: may want to ignore certain classes
		next if $ignore eq 'singleton' and $class eq $group->{$class}->[0];
		
		next if $ignore eq 'class'     and $class ne $group->{$class}->[0];
	    }
	    
	    #choose smallest class exceeding threshold and
	    #highest percent when same size
	    
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
		$bstclass = '*' #some non-gaps
	    }
	} else {
	    $bstclass = '*' #wildcard
	}
	#warn "DECIDE [$i] '$bstclass' $bstscore [$group->{$bstclass}->[0]]\n";
	$consensus .= $group->{$bstclass}->[0];
    }
    \$consensus;
}

sub new {
    my $type = shift;
    #warn "${type}::new() (@_)\n";
    if (@_ < 5) {
	die "${type}::new() missing arguments\n";
    }
    my ($from, $to, $tally, $group, $threshold, $ignore) = @_;

    if ($threshold < 50 or $threshold > 100) {
	die "${type}::new() threshold '$threshold\%' outside valid range [50..100]\n";
    }

    my $self = { %Bio::MView::Align::Sequence::Template };

    $self->{'id'}        = "consensus/$threshold\%";
    $self->{'type'}      = 'consensus';
    $self->{'from'}      = $from;
    $self->{'to'}        = $to;
    $self->{'threshold'} = $threshold;
    $self->{'group'}     = $group;

    my $string = consensus($tally, $group, $threshold, $ignore);

    #encode the new "sequence"
    $self->{'string'} = new Bio::MView::Sequence;
    $self->{'string'}->set_find_pad('\.');
    $self->{'string'}->set_find_gap('\.');
    $self->{'string'}->set_pad('.');
    $self->{'string'}->set_gap('.');
    $self->{'string'}->insert([$string, $from, $to]);
    
    bless $self, $type;

    $self->reset_display;

    $self;
}

sub color_by_type {
    my $self = shift;
    my %par = @_;

    #overkill!
    return unless $self->{'type'} eq 'consensus';

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Bio::MView::Align::Sequence::Default_Colormap
	unless defined $par{'colormap'};
    $par{'colormap2'}= $Default_Colormap
	unless defined $par{'colormap2'};

    my ($color, $end, $i, $cg, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    #warn "color_by_type($self) 1=$par{'colormap'} 2=$par{'colormap2'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$cg = $self->{'string'}->raw($i);
	
	#warn "[$i]= $cg\n";

	#white space: no color
	next    if $self->{'string'}->is_space($cg);

	#gap: gapcolour
	if ($self->{'string'}->is_non_char($cg)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#use symbol color/wildcard colour
	@tmp = $self->get_color_type($cg,
				     $par{'colormap'},
				     $par{'colormap2'});
	
	if (@tmp) {
	    if ($par{'css1'}) {
		push @$color, $i, 'class' => $tmp[1];
	    } else {
		push @$color, $i, 'color' => $tmp[0];
	    }
	} else {
	    push @$color, $i, 'color' => $par{'symcolor'};
	}
    }
    
    $self->{'display'}->{'paint'}  = 1;
    $self;
}

sub color_by_identity {
    my ($self, $othr) = (shift, shift);    #ignore second arg
    my %par = @_;

    #overkill!
    return unless $self->{'type'} eq 'consensus';

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Bio::MView::Align::Sequence::Default_Colormap
	unless defined $par{'colormap'};
    $par{'colormap2'}= $Default_Colormap
	unless defined $par{'colormap2'};

    my ($color, $end, $i, $cg, @tmp) = ($self->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $par{'symcolor'};
    
    #warn "color_by_identity($self, $othr) 1=$par{'colormap'} 2=$par{'colormap2'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$cg = $self->{'string'}->raw($i);

	#white space: no colour
	next    if $self->{'string'}->is_space($cg);
					 
	#gap: gapcolour
	if ($self->{'string'}->is_non_char($cg)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#consensus group symbol is singleton: choose colour
	if (exists $Group->{$self->{'group'}}->[2]->{$cg}) {
	    if (keys %{$Group->{$self->{'group'}}->[2]->{$cg}} == 1) {

		@tmp = $self->get_color_identity($cg, $par{'colormap'});

		if (@tmp) {
		    if ($par{'css1'}) {
			push @$color, $i, 'class' => $tmp[1];
		    } else {
			push @$color, $i, 'color' => $tmp[0];
		    }
		} else {
		    push @$color, $i, 'color' => $par{'symcolor'};
		}

		next;
	    }
	}
	
	#symbol not in consensus group: use contrast colour
	push @$color, $i, 'color' => $par{'symcolor'};
    }
    
    $self->{'display'}->{'paint'} = 1;
    $self;
}

#this is analogous to Bio::MView::Align::Row::Sequence::color_by_identity()
#but the roles of self (consensus) and other (sequence) are reversed.
sub color_by_consensus_sequence {
    my ($self, $othr) = (shift, shift);

    return unless $othr;
    return unless $othr->{'type'} eq 'sequence';

    my %par = @_;

    die "${self}::color_by_consensus_sequence() length mismatch\n"
	unless $self->length == $othr->length;

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Bio::MView::Align::Sequence::Default_Colormap
	unless defined $par{'colormap'};
    $par{'colormap2'}= $Default_Colormap
	unless defined $par{'colormap2'};

    my ($color, $end, $i, $cg, $cs, $c, @tmp) = ($othr->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    #warn "color_by_consensus_sequence($self, $othr) 1=$par{'colormap'} 2=$par{'colormap2'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$cg = $self->{'string'}->raw($i); $cs = $othr->{'string'}->raw($i);

	#warn "[$i]= $cg <=> $cs\n";

	#white space: no colour
	next    if $self->{'string'}->is_space($cs);
					 
	#gap: gapcolour
	if ($self->{'string'}->is_non_char($cs)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#symbols in consensus group are stored upcased
	$c = uc $cs;

	#symbol in consensus group: choose colour
	if (exists $Group->{$self->{'group'}}->[1]->{$c}) {
	    if (exists $Group->{$self->{'group'}}->[1]->{$c}->{$cg}) {

		#colour by sequence symbol
		@tmp = $self->get_color_consensus_sequence($cs, $cg,
							   $par{'colormap'},
							   $par{'colormap2'});

		if (@tmp) {
		    if ($par{'css1'}) {
			push @$color, $i, 'class' => $tmp[1];
		    } else {
			push @$color, $i, 'color' => $tmp[0];
		    }
		} else {
		    push @$color, $i, 'color' => $par{'symcolor'};
		}

		next;
	    }
	}

        #symbol not in consensus group: use contrast colour
	push @$color, $i, 'color' => $par{'symcolor'};
    }
    
    $othr->{'display'}->{'paint'} = 1;
    $self;
}


#this is analogous to Bio::MView::Align::Row::Sequence::color_by_identity()
#but the roles of self (consensus) and other (sequence) are reversed.
sub color_by_consensus_group {
    my ($self, $othr) = (shift, shift);

    return unless $othr;
    return unless $othr->{'type'} eq 'sequence';

    my %par = @_;

    die "${self}::color_by_consensus_group() length mismatch\n"
	unless $self->length == $othr->length;

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Bio::MView::Align::Sequence::Default_Colormap
	unless defined $par{'colormap'};
    $par{'colormap2'}= $Default_Colormap
	unless defined $par{'colormap2'};

    my ($color, $end, $i, $cg, $cs, $c, @tmp) = ($othr->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    #warn "color_by_consensus_group($self, $othr) 1=$par{'colormap'} 2=$par{'colormap2'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$cg = $self->{'string'}->raw($i); $cs = $othr->{'string'}->raw($i);

	#warn "[$i]= $cg <=> $cs\n";
	
	#no sequence symbol: whitespace: no colour
	next    if $self->{'string'}->is_space($cs);

	#gap or frameshift: gapcolour
	if ($self->{'string'}->is_non_char($cs)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#symbols in consensus group are stored upcased
	$c = uc $cs;

	#symbol in consensus group: choose colour
	if (exists $Group->{$self->{'group'}}->[1]->{$c}) {
	    if (exists $Group->{$self->{'group'}}->[1]->{$c}->{$cg}) {

		#colour by consensus group symbol
		#note: both symbols passed; colormaps swapped
		@tmp = $self->get_color_consensus_group($cs, $cg,
							$par{'colormap'},
							$par{'colormap2'});
		if (@tmp) {
		    if ($par{'css1'}) {
			push @$color, $i, 'class' => $tmp[1];
		    } else {
			push @$color, $i, 'color' => $tmp[0];
		    }
		} else {
		    push @$color, $i, 'color' => $par{'symcolor'};
		}

		next;
	    }
	}
	
	#symbol not in consensus group: use contrast colour
	push @$color, $i, 'color' => $par{'symcolor'};
    }
    
    $othr->{'display'}->{'paint'} = 1;
    $self;
}


###########################################################################
1;

__DATA__
#label      => symbol  { member list }

[P1]
#Protein consensus: conserved physicochemical classes, derived from
#the Venn diagrams of: Taylor W. R. (1986). The classification of amino acid
#conservation. J. Theor. Biol. 119:205-218.
*            =>  .
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

[D1]
#DNA consensus: conserved ring types
#Ambiguous base R is A or G
#Ambiguous base Y is C or T or U
*            =>  .
A            =>  A       { A }
G            =>  G       { G }
C            =>  C       { C }
T            =>  T       { T }
U            =>  U       { U }
purine       =>  r       { A, G,  R }
pyrimidine   =>  y       { C, T, U,  Y }

[CYS]
#Protein consensus: conserved cysteines
*            =>  .
C            =>  C       { C }


###########################################################################
