# Copyright (C) 2015 Nigel P. Brown
# $Id: Header.pm,v 1.1 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Display::Header;

use Bio::MView::Display::Any;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Display::Any);

sub new {
    my $type = shift;
    my $self = new Bio::MView::Display::Any(@_);
    bless $self, $type;
}

sub next {
    my $self = shift;
    my ($html, $bold, $col, $gap, $pad, $lap, $ruler) = @_;
    my ($string, $length, $i, $start, $pos) = ([]);
	 					  
    return 0    if $self->{'cursor'} > $self->{'length'};

    $length = $self->{'length'} - $self->{'cursor'} +1; #length remaining
    $col = $length    if $col > $length;                #consume smaller amount
    $start = $self->{'start'} + $self->{'cursor'} -1;   #current real position

    #warn "($self->{'cursor'}, $col, $length, ($self->{'start'},$self->{'stop'}))\n";

    for ($i = 0; $i < $col; $i++) {
	$pos = $start + $i;     #real data position
	push @$string, ' ';
    }

    $string= $self->strip_html_repeats($string)    if $html;

    if ($html) {
	unshift @$string, '<STRONG>'           if $bold;
	unshift @$string, $self->{'prefix'}    if $self->{'prefix'};
	push    @$string, $self->{'suffix'}    if $self->{'suffix'};
	push    @$string, '</STRONG>'          if $bold;
    }

    $self->{'cursor'} += $col;

    return [ $start, $pos, join('', @$string) ];
}


###########################################################################
1;
