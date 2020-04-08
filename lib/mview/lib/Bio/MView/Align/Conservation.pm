# Copyright (C) 2015-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Align::Conservation;

use Bio::MView::Sequence::Base;
use Bio::MView::Align::Sequence;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Sequence);

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    die "${type}::new: missing arguments\n"  if @_ < 3;
    my ($from, $to, $string) = @_;

    #encode the new "sequence"
    my $sob = new Bio::MView::Sequence;
    $sob->set_find_pad(' '); $sob->set_pad(' ');
    $sob->set_find_gap(' '); $sob->set_gap(' ');
    $sob->insert([$string, $from, $to]);

    my $self = new Bio::MView::Align::Sequence('clustal', $sob);
    bless $self, $type;

    $self;
}

######################################################################
# public methods
######################################################################
#override
sub is_sequence { 0 }

###########################################################################
1;
