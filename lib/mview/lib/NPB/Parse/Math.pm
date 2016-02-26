# -*- perl -*-
# Copyright (C) 1996-2006 Nigel P. Brown
# $Id: Math.pm,v 1.2 2005/12/12 20:42:48 brown Exp $

###########################################################################
package NPB::Parse::Math;

use Exporter;

@ISA = Exporter;

@EXPORT = qw(
	     min
	     max
	     );

use strict;

###########################################################################
#arithmetic min() function
sub min {
    my ($a, $b) = @_;
    $a < $b ? $a : $b;
}

#arithmetic max() function
sub max {
    my ($a, $b) = @_;
    $a > $b ? $a : $b;
}

###########################################################################
1;
