# Copyright (C) 2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Display::Display;

use Bio::Util::Math qw(max);
use Bio::Util::System qw(vmstat);
use Bio::MView::Display::Panel;
use Bio::MView::Display::Out::Text;
use Bio::MView::Display::Out::HTML;

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";

    my $self = {};
    bless $self, $type;

    $self->{'par'}   = undef;
    $self->{'panel'} = [];

    $self->initialise_parameters(@_);
    $self->initialise_fieldwidths;

    $self;
}

######################################################################
# public methods
######################################################################
sub new_panel {
    my $self = shift;

    my $par = $self->{'par'};

    my $panel = new Bio::MView::Display::Panel($par, @_);
    push @{$self->{'panel'}}, $panel;

    return $panel;
}

sub render {
    my ($self, $mode) = (@_, undef);

    return  unless @{$self->{'panel'}};

    my $par = $self->{'par'};

    my ($posnwidths, $labelwidths) = $self->update_multipanel_fieldwidths;

    $par->{'dev'} = select_output($par, $mode, $self->{'par'}->{'stream'});

    while (@{$self->{'panel'}}) {  #consume panels

        my $pan = shift @{$self->{'panel'}};

        $par->{'dev'}->render_table_begin($par);
        $pan->render_panel($par, $posnwidths, $labelwidths);
        $par->{'dev'}->render_table_end;
    }
}

######################################################################
# private methods
######################################################################
sub initialise_parameters {
    my $self = shift;
    my $par = {@_};

    $par->{'width'} = 0
        unless exists $par->{'width'} and defined $par->{'width'};

    $par->{'bold'} = 0
        unless exists $par->{'bold'} and defined $par->{'bold'};

    $self->{'par'} = $par;
}

sub initialise_fieldwidths {
    my $self = shift;

    my $par = $self->{'par'};

    #initialise left/right posnwidths
    $self->{'posnwidths'}  = [0, 0];

    #initialise labelwidths
    $self->{'labelwidths'} = [];
    for (my $i=0; $i < @{$par->{'labelflags'}}; $i++) {
        $self->{'labelwidths'}->[$i] = 0;
    }
}

sub update_multipanel_fieldwidths {
    my $self = shift;

    return (undef, undef)  unless $self->{'par'}->{'register'};

    my $posnwidths  = $self->{'posnwidths'};
    my $labelwidths = $self->{'labelwidths'};

    #consolidate widths across multiple Panel objects
    foreach my $pan (@{$self->{'panel'}}) {

        #numeric left/right position widths
        for (my $i=0; $i < @$posnwidths; $i++) {
            $posnwidths->[$i] = max($posnwidths->[$i], $pan->posnwidth($i));
        }

        #labelwidths
        for (my $i=0; $i < @$labelwidths; $i++) {
            $labelwidths->[$i] = max($labelwidths->[$i], $pan->labelwidth($i));
        }
    }

    $self->{'posnwidths'}  = $posnwidths;
    $self->{'labelwidths'} = $labelwidths;

    return ($posnwidths, $labelwidths);
}

sub select_output {
    my ($par, $mode, $stm) = @_;

    if (! defined $mode) {
        $mode = 'html'  if $par->{'html'};
        $mode = 'text'  if !$par->{'html'};
    }

    return new Bio::MView::Display::Out::Text($stm)  if $mode eq 'text';
    return new Bio::MView::Display::Out::HTML($stm)  if $mode eq 'html';

    die "Panel: unknown output mode '§mode'\n";
}

###########################################################################
1;
