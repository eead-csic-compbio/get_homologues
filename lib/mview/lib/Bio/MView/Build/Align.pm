# Copyright (C) 1997-2006 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Build::Align;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Base;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Base);

######################################################################
# protected methods
######################################################################
#override: do all rows have same range?
sub test_if_aligned {
    my ($self, $lo, $hi) = @_;
    for (my $i=1; $i < @{$self->{'index2row'}}; $i++) { #from 1
        my $row = $self->{'index2row'}->[$i];
        my ($lo2, $hi2) = $self->get_range($row);
        #warn "$row ($lo2, $hi2)\n";
        return 0  if $lo != $lo2 or $hi != $hi2;
    }
    return 1;
}

#override
sub use_row {
    my ($self, $num, $nid, $sid) = @_;
    my $pat;

    #warn "use_row($num, $nid, $sid)\n";

    #first, check explicit keeplist and reference row
    foreach $pat (@{$self->{'keeplist'}}, $PAR->get('ref_id')) {

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
    foreach $pat (@{$self->{'skiplist'}}, $PAR->get('ref_id')) {

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

#override
sub map_id {
    my ($self, $id) = @_;
    $id = 1  if $id eq '0' or $id =~ /query/i;
    return $self->SUPER::map_id($id);
}

###########################################################################
1;
