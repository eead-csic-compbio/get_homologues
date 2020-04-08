# Copyright (C) 1996-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Util::Regexp;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);

###########################################################################
#integers
###########################################################################

use vars qw($RX_Uint $RX_Sint);

push @EXPORT, qw($RX_Uint $RX_Sint);

#unsigned
$RX_Uint  = '\+?\d+';

#signed
$RX_Sint  = '[+-]?\d+';

###########################################################################
#floats
###########################################################################

use vars qw($RX_Ureal $RX_Sreal);

push @EXPORT, qw($RX_Ureal $RX_Sreal);

#unsigned
$RX_Ureal = '\+?(?:\d+\.\d+|\d+\.|\d+|\.\d+)?(?:[eE][+-]?\d+)?';

#signed
$RX_Sreal = '[+-]?(?:\d+\.\d+|\d+\.|\d+|\.\d+)?(?:[eE][+-]?\d+)?';

###########################################################################
1;
