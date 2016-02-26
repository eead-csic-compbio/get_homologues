# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Any.pm,v 1.12 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Display::Any;

use Bio::MView::Display;
use Bio::MView::Display::Ruler;
use Bio::MView::Display::Header;

use strict;
use vars qw(%Template);

%Template =
    (
     'parent'  	=> undef, #parent sequence object
     'length'  	=> 0,     #length of parent sequence
     'start'   	=> 0,     #start position of parent sequence
     'stop'    	=> 0,     #stop position of parent sequence
     'string'  	=> undef, #alternative sequence to parent
     'type'    	=> '',    #kind of aligned object
     'label0'  	=> '',    #first label string
     'label1'  	=> '',    #second label string
     'label2'  	=> '',    #third label string
     'label3'  	=> '',    #fourth label string
     'label4'  	=> '',    #fifth label string
     'label5'  	=> '',    #sixth label string
     'label6'  	=> '',    #seventh label string
     'number'  	=> 0,     #marginal numbers flag
     'url'     	=> '',    #URL to associate with label1
     'cursor'  	=> 0,     #current position in display stream, 1-based
     'paint'   	=> 0,     #overlap behaviour for this row
     'prefix' 	=> '',    #prefix whole row with string
     'suffix' 	=> '',    #suffix whole row with string

     'r_map'  	=> undef,
     'r_color' 	=> undef,
     'r_class' 	=> undef,
     'r_prefix'	=> undef,
     'r_suffix'	=> undef,
     'r_sym'   	=> undef,
     'r_url'   	=> undef,
     'r_case'  	=> undef,
    );

sub new {
    my $type = shift;
    my ($parent, $hash) = @_;
    my ($range, $i, $j);
    my $self = { %Template };
    bless $self, $type;

    $self->{'parent'}  	= $parent;
    $self->{'length'}  	= $parent->length;
    $self->{'type'}    	= $hash->{'type'};
    $self->{'string'}  	= $hash->{'sequence'} if exists $hash->{'sequence'};
    $self->{'label0'}  	= $hash->{'label0'}   if exists $hash->{'label0'};
    $self->{'label1'}  	= $hash->{'label1'}   if exists $hash->{'label1'};
    $self->{'label2'}  	= $hash->{'label2'}   if exists $hash->{'label2'};
    $self->{'label3'}  	= $hash->{'label3'}   if exists $hash->{'label3'};
    $self->{'label4'}  	= $hash->{'label4'}   if exists $hash->{'label4'};
    $self->{'label5'}  	= $hash->{'label5'}   if exists $hash->{'label5'};
    $self->{'label6'}  	= $hash->{'label6'}   if exists $hash->{'label6'};
    $self->{'number'}  	= $hash->{'number'}   if exists $hash->{'number'};
    $self->{'url'}     	= $hash->{'url'}      if exists $hash->{'url'};
    $self->{'paint'}   	= $hash->{'paint'}    if exists $hash->{'paint'};
    $self->{'prefix'} 	= $hash->{'prefix'}   if exists $hash->{'prefix'};
    $self->{'suffix'} 	= $hash->{'suffix'}   if exists $hash->{'suffix'};
 
    $self->{'r_map'}   	= [];
    $self->{'r_color'} 	= [];
    $self->{'r_class'} 	= [];
    $self->{'r_prefix'}	= [];
    $self->{'r_suffix'}	= [];
    $self->{'r_sym'}   	= [];
    $self->{'r_url'}   	= [];
    $self->{'r_case'}  	= [];

    if (defined $self->{'string'}) {
	my $seqlen = $self->{'string'}->length;

	#are 'parent' and 'string' the same length?
	if ($seqlen != $self->{'length'}) {
            #warn "parent:   [@{[$self->{'parent'}->{'string'}->string]}]\n";
            #warn "sequence: [@{[$self->{'string'}->string]}]\n";
	    die "$type new() parent/sequence length mismatch ($self->{'length'}, $seqlen)\n";
	}

	# no idea what this is/was for... 1/10/98
	#starting numbering wrt 'string' or 'parent'?
	if (exists $hash->{'start'}) {
	    $self->{'start'} = $hash->{'start'};
	    $self->{'stop'}  = $hash->{'stop'};
	} else {
	    $self->{'start'} = $parent->{'start'};
	    $self->{'stop'}  = $parent->{'stop'};
	}

    } else {
	#starting numbering wrt 'parent'
	$self->{'start'} = $parent->{'start'};
        $self->{'stop'}  = $parent->{'stop'};
    }

    #ensure 'range' attribute has a value
    $hash->{'range'} = []  unless exists $hash->{'range'};

    #populate the range map
    foreach $range (@{$self->parse_range($hash->{'range'})}) {
	for ($i = $range->{'pos'}->[0];
	     $i <= $range->{'pos'}->[$#{$range->{'pos'}}];
	     $i++) {
	    foreach $j (keys %$range) {
		next if $j eq 'pos';
		$self->{"r_$j"}->[$i] = $range->{$j};
	    }
            $self->{'r_map'}->[$i]++;    #count hits per position
	}
    }

    #Universal::vmstat("Display::Any::new done");
    $self;
}

sub parse_range {
    my ($self, $data) = @_;
    my ($state, @list, $attr, $i) = (2);
    
    for ($i=0; $i < @$data; $i++) {

	#position data
	if ($data->[$i] =~ /^\d+$/) {
	    if ($state == 2) {
		#warn "SAVE   [$i] $data->[$i]\n";
		#begin new range: save old attribute hash, if useful
		if (defined $attr and keys %$attr > 1) {
		    if ($attr->{'pos'}->[0] > $attr->{'pos'}->[1]) {
			warn "bad attribute range '@{$attr->{'pos'}}'\n";
		    } else {
			push @list, $attr;
		    }
		}
		$attr  = {};   #reset attribute hash
		$state = 0;    #ready for first pos
	    }
	    if ($state == 0) {
		#warn "0->1   [$i] $data->[$i]\n";
		#first position
		$attr->{'pos'}->[0] = $data->[$i];
		$state = 1;    #ready for second pos, if any
		next;
	    }
	    if ($state == 1) {
		#warn "1->2   [$i] $data->[$i]\n";
		#second position
		$attr->{'pos'}->[1] = $data->[$i];
		$state = 2;
		next;
	    }
	    die "$self new() never get here!\n";
	}

	#first position wasn't given: bad attribute list
	if ($state == 0) {
	    warn "$self new() bad attribute list '@$data'\n";
	    last;
	}
	
	#second position wasn't given: make one and continue
	if ($state == 1) {
	    $attr->{'pos'}->[1] = $attr->{'pos'}->[0];
	    $state = 2;
	}
	
	#position range attributes: key=>value pair
	if (grep /^$data->[$i]$/, qw(color class prefix suffix sym url case)) {
	    #warn "ATTR   [$i] $data->[$i] $data->[$i+1]\n";
	    $attr->{$data->[$i]} = $data->[$i+1]  if defined $data->[$i+1];
	    $i++;    #jump value
	    next;
	}
    
	warn "$self new() unknown range attribute '$data->[$i]'\n";
	$i++;    #skip likely value
    }

    #save last iteration
    if ($state == 2 and defined $attr) {
	#warn "SAVE   LAST\n";
	#save old attribute hash, if useful
	if (defined $attr and keys %$attr > 1) {
	    if ($attr->{'pos'}->[0] > $attr->{'pos'}->[1]) {
		warn "bad attribute range '@{$attr->{'pos'}}'\n";
	    } else {
		push @list, $attr;
	    }
	}
    }

    #resort ranges by (1) lowest 'start', and (2) highest 'stop'
    return [
	    sort {
		my $c = $a->{'pos'}->[0] <=> $b->{'pos'}->[0];
		return  $c if $c != 0;
		return  $b->{'pos'}->[1] <=> $a->{'pos'}->[1];
	    } @list
	   ];
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

#force break of back reference to allow garbage collection
sub free {
    local $_;
    #warn "FREE $_[0]\n";
    $_[0]->{'parent'} = undef;
}

sub type   { $_[0]->{'type'} }
sub label0 { $_[0]->{'label0'} }
sub label1 { $_[0]->{'label1'} }
sub label2 { $_[0]->{'label2'} }
sub label3 { $_[0]->{'label3'} }
sub label4 { $_[0]->{'label4'} }
sub label5 { $_[0]->{'label5'} }
sub label6 { $_[0]->{'label6'} }

sub col { 
    defined $_[0]->{'string'} ?
	$_[0]->{'string'}->col($_[1]) :
	    $_[0]->{'parent'}->{'string'}->col($_[1]);
}

sub print {
    my $self = shift;
    local $_;
    foreach (sort keys %$self) {
	printf "%15s => %s\n", $_, $self->{$_};
    }
    $self;
}

#forward/reverse oriented rulers/sequences
sub reset { $_[0]->{'cursor'} = 1 }

sub done {
    return 1    if $_[0]->{'cursor'} > $_[0]->{'length'};
    return 0;
}

sub next {
    my $self = shift;
    #warn "${self}::next(@_)\n";
    my ($mode, $html, $bold, $col, $gap, $pad, $lap, $ruler) = @_;
    my (@string, $length, $i, $start, $pos, $c) = ();
	 					  
    return 0    if $self->{'cursor'} > $self->{'length'};

    $length = $self->{'length'} - $self->{'cursor'} +1; #length remaining
    $col = $length    if $col > $length;                #consume smaller amount
    $start = $self->{'start'} + $self->{'cursor'} -1;   #current real position

    $col++;

    #warn "($self->{'cursor'}, $col, $length, ($self->{'start'},$self->{'stop'}))\n";

    for ($i = 1; $i < $col; $i++, $self->{'cursor'}++) {
	
	$pos = $start + $i;     #real data position
	
	#warn "$self->{'cursor'}, $pos ", mapcount=$self->{'r_map'}, "\n";

	#any range specifications?
	if (defined $self->{'r_map'}->[$self->{'cursor'}]) {

	    #overlapped position if 'paint' mode is off
	    if (! $self->{'paint'} and
		$self->{'r_map'}->[$self->{'cursor'}] > 1) {

		if ($html and length $lap > 1) {
		    push @string, "<SPAN style=\"color:$lap\">";
		    push @string, $self->col($self->{'cursor'});
		    push @string, "</SPAN>";
		} else {
		    push @string, $lap;
		}
		next;
	    }

	    #non-overlapped position, or 'paint' mode is on
	    if (defined $self->{'r_sym'}->[$self->{'cursor'}]) {
		#caller specified symbol
		$c = $self->{'r_sym'}->[$self->{'cursor'}];
	    } elsif (defined $self->{'r_case'}->[$self->{'cursor'}]) {
		#caller specified case
		$c = $self->{'r_case'}->[$self->{'cursor'}];
		$c = uc $self->col($self->{'cursor'})  if $c eq 'uc';
		$c = lc $self->col($self->{'cursor'})  if $c eq 'lc';
	    } else {
		#default to sequence symbol
		$c = $self->col($self->{'cursor'});
	    }
	    push @string, $self->html_wrap($self->{'cursor'}, $c)  if $html;

	} elsif ($mode eq 'gap') {
	    #class Subrange: use gap character
	    push @string, $gap;
	} else { #mode eq 'nogap'
	    #class Sequence: use sequence character
	    push @string, $self->col($self->{'cursor'});
        }
    }
    
    @string = @{ $self->strip_html_repeats(\@string) }    if $html;

    if ($html) {
	unshift @string, '<STRONG>'         if $bold;
	unshift @string, $self->{'prefix'}  if $self->{'prefix'};
	push    @string, $self->{'suffix'}  if $self->{'suffix'};
	push    @string, '</STRONG>'        if $bold;
    }

    return [ $start, $pos, join('', @string) ];
}

sub html_wrap {
    my ($self, $i, $c) = @_;

    my @string = ();

    #open new url?
    push @string, "<A HREF=\"$self->{'r_url'}->[$i]\">"
	if defined $self->{'r_url'}->[$i];

    #change color?
    push @string, "<SPAN style=\"color:$self->{'r_color'}->[$i]\">"
	if defined $self->{'r_color'}->[$i] and
	    ! defined $self->{'r_class'}->[$i];
    
    #css1 class?
    push @string, "<SPAN CLASS=$self->{'r_class'}->[$i]>"
	if defined $self->{'r_class'}->[$i];
    
    #prefix
    push @string, $self->{'r_prefix'}->[$i]
        if defined $self->{'r_prefix'}->[$i];
    
    #embedded character
    push @string, $c;

    #suffix
    push @string, $self->{'r_suffix'}->[$i]
        if defined $self->{'r_suffix'}->[$i];
    
    #unchange css1 class?
    push @string, "</SPAN>"  if defined $self->{'r_class'}->[$i];
    
    #unchange color?
    push @string, "</SPAN>"
	if defined $self->{'r_color'}->[$i] and
	    ! defined $self->{'r_class'}->[$i];
    
    #close URL?
    push @string, "</A>"     if defined $self->{'r_url'}->[$i];
    
    @string;
}

sub strip_html_repeats {
    my ($self, $list) = @_;
    my ($new, $i, $limit) = [];

    return $list  if @$list < 4;

    #warn " IN=[", @$list, "]\n";

    #use 1 element lookahead
    $limit = @$list;

    for ($i=0; $i<$limit; $i++) {
	
	if ($list->[$i] =~ /^</) {
	    #looking at: <tag> or </tag>

	    #empty stack
	    if (@$new < 1) {
		#warn "NEW(@{[scalar @$new]}) lex[$i] <-- $list->[$i]\n";
		push @$new, $list->[$i];
		#warn "[", @$new, "]\n";
		next;
	    }
	    
	    #non-empty stack: different tag
	    if ($list->[$i] ne $new->[$#$new]) {
		#warn "NEW(@{[scalar @$new]}) lex[$i] <-- $list->[$i]\n";
		push @$new, $list->[$i];
		#warn "[", @$new, "]\n";
		next;
	    }

	    #non-empty stack: same tag
	    #warn "NOP(@{[scalar @$new]}) lex[$i]     $new->[$#$new]\n";
	    #warn "[", @$new, "]\n";
	    
	} else {
	    #non-tag
	    #warn "ARG(@{[scalar @$new]})    [$i] <-- $list->[$i]\n";
	    push @$new, $list->[$i];
	    #warn "[", @$new, "]\n";
	}
    }
    
    #warn "SR1=[", @$new, "]\n";

    $new = $self->strip_html_tandem_repeats($new);

    #warn "SR2=[", @$new, "]\n\n";

    $new;
}

sub strip_html_tandem_repeats {
    my ($self, $list) = @_;
    my ($new, @mem, $i, $limit) = ([]);

    return $list  if @$list < 6;

    #warn " IN=[", @$list, "]\n";

    #use 1 element lookahead
    $limit = @$list;

    for ($i=0; $i<$limit; $i++) {

        #looking at: tag
        if ($list->[$i] =~ /^<[^\/]/) {

            #empty stack
            if (@mem < 1) {
		#warn "NEW(@{[scalar @mem]}) lex[$i] <-- $list->[$i]\n";
                push @mem, $list->[$i];
                #scan until next closing tag
                $i = _close_tag($list, $limit, ++$i, \@mem);
                next;
            }

            #non-empty stack: different tag
            if ($list->[$i] ne $mem[0]) {
		#warn "NEW(@{[scalar @mem]}) lex[$i] <-- $list->[$i]\n";
                push @$new, @mem;
                @mem = ();
                redo;
            }

            #non-empty stack: same tag
	    #warn "POP(@{[scalar @mem]}) lex[$i] --> $mem[$#mem]\n";
            pop @mem;
            #scan until next closing tag
            $i = _close_tag($list, $limit, ++$i, \@mem);
            next;

        } else {
            #non-tag
	    #warn "ARG(@{[scalar @mem]}) lex[$i] <-- $list->[$i]\n";
            push @$new, @mem, $list->[$i];
            @mem = ();
        }
    }

    push @$new, @mem    if @mem;

    #warn "SR2=[", @$new, "]\n\n";

    $new;
}

sub _close_tag {
    my ($list, $limit, $i, $mem) = @_;
    my $tag = 0;
    while ($i < $limit) {
	
	if ($list->[$i] =~ /^<[^\/]/) {
	    #inner tag
	    $tag++;
	    #warn "I  tag[$i] <-- $list->[$i]\n";
	    push @$mem, $list->[$i++];
	    next;
	}
	
	if ($list->[$i] =~ /^<\//) {
	    if ($tag) {
		#inner /tag
		$tag--;
		#warn "I /tag[$i] <-- $list->[$i]\n";
		push @$mem, $list->[$i++];
		next;
	    }
	    #outer /tag
	    #warn "O /tag[$i] <-- $list->[$i]\n";
	    push @$mem, $list->[$i];
	    last;
	}
	
	#datum
	#warn "  data[$i] <-- $list->[$i]\n";
	push @$mem, $list->[$i++];
    }
    $i;
}


###########################################################################
package Bio::MView::Display::Sequence;

@Bio::MView::Display::Sequence::ISA = qw(Bio::MView::Display::Any);

sub next { my $self = shift; $self->SUPER::next('nogap', @_) }


###########################################################################
package Bio::MView::Display::Subrange;

@Bio::MView::Display::Subrange::ISA = qw(Bio::MView::Display::Any);

sub next { my $self = shift; $self->SUPER::next('gap', @_) }


###########################################################################
1;
