# Copyright (C) 1997-2006 Nigel P. Brown
# $Id: Align.pm,v 1.12 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Build::Align;

use Bio::MView::Build;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build);

sub use_row {
    my ($self, $num, $nid, $sid) = @_;
    my $pat;

    #warn "use_row($num, $nid, $sid)\n";

    #first, check explicit keeplist and reference row
    foreach $pat (@{$self->{'keeplist'}}, $self->{'ref_id'}) {

	#Align subclass only
	return 1  if $pat eq '0'     and $num == 1;
	return 1  if $pat eq 'query' and $num == 1;

	#look at row number
	return 1  if $nid eq $pat;      #major
	
	#look at identifier
	return 1  if $sid eq $pat;      #exact match
	if ($pat =~ /^\/(.*)\/$/) {     #regex match (case insensitive)
	    return 1  if $sid =~ /$1/i;
	}
    }

    #second, check skiplist and reference row
    foreach $pat (@{$self->{'skiplist'}}, $self->{'ref_id'}) {

	#Align subclass only
	return 0  if $pat eq '0'     and $num == 1;
	return 0  if $pat eq 'query' and $num == 1;

	#look at row number
	return 0  if $nid eq $pat;      #major
	
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
    $ref = 1  if $ref eq '0' or $ref =~ /query/i;
    $self->SUPER::map_id($ref);
}


###########################################################################
1;
