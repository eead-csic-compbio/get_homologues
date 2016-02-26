# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Sequence.pm,v 1.26 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Align::Sequence;

use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Align::Row;

use strict;
use vars qw(@ISA
	    $Default_PRO_Colormap $Default_DNA_Colormap
            $Default_FIND_Colormap $Default_Colormap
	    %Template);

@ISA = qw(Bio::MView::Align::Row);

$Default_PRO_Colormap  = 'P1';    #default protein colormap name
$Default_DNA_Colormap  = 'D1';    #default nucleotide colormap name
$Default_FIND_Colormap = 'FIND';  #default find pattern colormap name
$Default_Colormap = $Default_PRO_Colormap;

%Template = 
    (
     'id'      => undef,     #identifier
     'type'    => undef,     #information about own subtype
     'from'    => undef,     #start number of sequence
     'string'  => undef,     #alignment string
     'display' => undef,     #hash of display parameters
    );

sub new {
    my $type = shift;
    #warn "${type}::new() (@_)\n";
    if (@_ < 2) {
	die "${type}::new() missing arguments\n";
    }
    my ($id, $string, $subtype) = (@_, 'sequence');

    my $self = { %Template };

    $self->{'id'}     = $id;
    $self->{'type'}   = $subtype;
    $self->{'from'}   = $string->lo;
    $self->{'string'} = $string;

    bless $self, $type;

    $self->reset_display;

    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

sub print {
    sub _format {
	my ($self, $k, $v) = @_;
	$v = 'undef' unless defined $v;
	$v = "'$v'" if $v =~ /^\s*$/;
	return sprintf("  %-15s => %s\n", $k, $v)
    }
    my $self = shift;
    warn "$self\n";
    map { warn $self->_format($_, $self->{$_}) } sort keys %{$self};
    $self;
}

sub id       { $_[0]->{'id'} }
sub seqobj   { $_[0]->{'string'} }
sub string   { $_[0]->{'string'}->string }
sub sequence { $_[0]->{'string'}->sequence }
sub from     { $_[0]->{'from'} }
sub length   { $_[0]->{'string'}->length }
sub seqlen   { $_[0]->{'string'}->seqlen }

sub reset_display {
    $_[0]->{'display'} =
	{
	 'type'     => 'sequence',
	 'label1'   => $_[0]->{'id'},
	 'sequence' => $_[0]->{'string'},
	 'range'    => [],
	};
    $_[0];
}

sub get_color {
    my ($self, $c, $map) = @_;
    my ($index, $color, $trans);

    #warn "get_color: $c, $map";

    #set transparent(T)/solid(S)
    if (exists $Bio::MView::Align::Colormaps->{$map}->{$c}) {
	$trans = $Bio::MView::Align::Colormaps->{$map}->{$c}->[1];
	$index = $Bio::MView::Align::Colormaps->{$map}->{$c}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $map\{$c} [$index] [$color] [$trans]\n";
	
	return ($color, "$trans$index");
    }

    #wildcard colour
    if (exists $Bio::MView::Align::Colormaps->{$map}->{'*'}) {

	$trans = $Bio::MView::Align::Colormaps->{$map}->{'*'}->[1];
	$index = $Bio::MView::Align::Colormaps->{$map}->{'*'}->[0];
	$color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $map\{'*'} [$index] [$color] [$trans]\n";

	return ($color, "$trans$index");
    }

    #preset colour name in $map, used for string searches
    #where all matches should be same colour
    if (exists $Bio::MView::Align::Palette->[0]->{$map}) {

	$trans = 'S';
	$index = $Bio::MView::Align::Palette->[0]->{$map};
        $color = $Bio::MView::Align::Palette->[1]->[$index];

	#warn "$c $map\{$c} [$index] [$color] [$trans]\n";

	return ($color, "$trans$index");
    }

    return 0;    #no match
}

sub color_special {
    my $self = shift;
    my %par = @_;

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};

    #locate a 'special' colormap'
    my ($size, $map) = (0);
    foreach $map (keys %$Bio::MView::Align::Colormaps) {
	if ($self->{'id'} =~ /$map/i) {
	    if (length($&) > $size) {
		$par{'colormap'} = $map;
		$size = length($&);
	    }
	}
    }
    return unless $size;

    #warn $par{'colormap'};

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$c = $self->{'string'}->raw($i);
	
	#warn "$self->{'id'}  [$i]= $c\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c);

	#gap: gapcolour
	if ($self->{'string'}->is_non_char($c)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#use symbol color/wildcard colour
	@tmp = $self->get_color($c, $par{'colormap'});

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

sub color_by_find_block {
    my $self = shift;
    my %par = @_;

    return  unless $self->{'type'} eq 'sequence';

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};
    $par{'find'}     = ''
        unless defined $par{'find'};

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    my %mindex = ();
    foreach my $block (@{$self->{string}->findall($par{patterns},
						  $par{mapsize})}) {
	$mindex{$block->[1]} = $block->[0];
    }

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {
        
	$c = $self->{'string'}->raw($i);
	
	#warn "[$i]= $c\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c);
        
	#gap or frameshift: gapcolour
	if ($self->{'string'}->is_non_char($c)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
        if (exists $mindex{$i}) {
            #use symbol color/wildcard colour
            @tmp = $self->get_color($mindex{$i}, $par{'colormap'});
        } else {
            @tmp = ();
        }
	
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
    
    $self->{'display'}->{'paint'} = 1;
    $self;
}

sub color_by_type {
    my $self = shift;
    my %par = @_;

    return  unless $self->{'type'} eq 'sequence';

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$c = $self->{'string'}->raw($i);
	
	#warn "[$i]= $c\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c);

	#gap or frameshift: gapcolour
	if ($self->{'string'}->is_non_char($c)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#use symbol color/wildcard colour
	@tmp = $self->get_color($c, $par{'colormap'});
	
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
    my ($self, $othr) = (shift, shift);
    my %par = @_;

    return  unless $self->{'type'} eq 'sequence';
    return  unless defined $othr;

    die "${self}::color_by_identity() length mismatch\n"
	unless $self->length == $othr->length;

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Default_Colormap
	unless defined $par{'colormap'};

    my ($color, $end, $i, $c1, $c2, @tmp) = ($self->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$c1 = $self->{'string'}->raw($i); $c2 = $othr->{'string'}->raw($i);

	#warn "[$i]= $c1 <=> $c2\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c1);
					 
	#gap of frameshift: gapcolour
	if ($self->{'string'}->is_non_char($c1)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}

	#same symbol when upcased: use symbol/wildcard color
	if (uc $c1 eq uc $c2) {

            @tmp = $self->get_color($c1, $par{'colormap'});

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

	#different symbol: use contrast colour
	#push @$color, $i, 'color' => $par{'symcolor'};
	
	#different symbol: use wildcard colour
	@tmp = $self->get_color('*', $par{'colormap'});
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

sub set_identity {
    #warn "Bio::MView::Align::Sequence::set_identity(@_)\n";
    my ($self, $ref, $mode) = @_;
    my $val = $self->compute_identity_to($ref, $mode);
    $self->set_display('label4'=>sprintf("%.1f%%", $val));
}

#compute percent identity to input reference object. Normalisation
#depends on the mode argument: 'reference' divides by the reference
#sequence length,'aligned' by the aligned region length (like blast),
#and 'hit' by the hit sequence. The last is the same as 'aligned' for
#blast, but different for multiple alignments like clustal.
sub compute_identity_to {
    #warn "Bio::MView::Align::Sequence::compute_identity_to(@_)\n";
    my ($self, $othr, $mode) = @_;

    return 0      unless defined $othr;
    return 100.0  if $self == $othr;  #always 100% identical to self

    die "${self}::compute_identity_to() length mismatch\n"
	unless $self->length == $othr->length;
    
    my ($sum, $len) = (0, 0);
    my $end = $self->length +1;

    for (my $i=1; $i<$end; $i++) {
	my $cnt = 0;
	
	my $c1 = $self->{'string'}->raw($i);
	my $c2 = $othr->{'string'}->raw($i);

	#at least one must be a sequence character
	$cnt++  if $self->{'string'}->is_char($c1);
	$cnt++  if $self->{'string'}->is_char($c2);
	next  if $cnt < 1;

        #standardize case
        $c1 = uc $c1; $c2 = uc $c2;

	#ignore terminal gaps in the *first* sequence
	$len++  unless $self->{'string'}->is_terminal_gap($c1);

        #ignore unknown character: contributes to length only
        next  if $c1 eq 'X' or $c2 eq 'X';

	$sum++  if $c1 eq $c2;
	#warn "[$i] $c1 $c2 : $cnt => $sum / $len\n";
    }
    
    #normalise identities
    my $norm = 0;
    if ($mode eq 'reference') {
	$norm = $othr->seqlen;
	#warn "ref norm= $norm\n";
    } elsif ($mode eq 'aligned') {
	$norm = $len;
	#warn "aln norm= $norm\n";
    } elsif ($mode eq 'hit') {
	$norm = $self->seqlen;
	#warn "hit norm= $norm\n";
    }
    #warn "identity $self->{'id'} = $sum/$norm\n";

    return ($sum = 100 * ($sum + 0.0) / $norm)    if $norm > 0;
    return 0;
}


###########################################################################
1;
