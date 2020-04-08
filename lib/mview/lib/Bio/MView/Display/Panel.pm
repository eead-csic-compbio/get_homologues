# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Display::Panel;

use Bio::Util::Math qw(max);
use Bio::Util::System qw(vmstat);
use Bio::MView::Display::Sequence;
use Bio::MView::Display::Ruler;
use Bio::MView::Display::Out::Text qw($LJUST $RJUST);

my $MINPOSWIDTH = 2;  #smallest ruler left/right position width
my $HSPACE      = 1;  #extra spaces between columns

my %Known_Track_Types = (
    'ruler'    => 1,
    'sequence' => 1,
    );

sub new {
    my $type = shift;
    die "${type}::new(@_): missing arguments\n"  if @_ < 6;
    my ($par, $headers, $forwards, $length, $start, $stop) = @_;

    my $self = {};
    bless $self, $type;

    $self->{'par'}      = $par;       #parameters
    $self->{'header'}   = $headers;   #list of header strings
    $self->{'forwards'} = $forwards;  #alignment orientation
    $self->{'length'}   = $length;    #alignment length
    $self->{'start'}    = $start;     #alignment start position
    $self->{'stop'}     = $stop;      #alignment stop position

    $self->{'track'}    = [];  #display objects

    #warn "${type}::new: start= $start, stop= $stop, (fw= $forwards)";

    $self->initialise_fieldwidths;

    $self;
}

######################################################################
# public methods
######################################################################
sub length     { return $_[0]->{'length'} }
sub forwards   { return $_[0]->{'forwards'} }
sub posnwidth  { return $_[0]->{'posnwidths'}->[$_[1]] }
sub labelwidth { return $_[0]->{'labelwidths'}->[$_[1]] }

sub append {
    my $self = shift;

    #handle data rows
    foreach my $row (@_) {

        unless (exists $row->{'type'}) {
            warn "${self}::append: missing type in '$row'\n";
            next;
        }

        my $type = $row->{'type'};

        unless (exists $Known_Track_Types{$type}) {
            warn "${self}::append: unknown alignment type '$type'\n";
            next;
        }

        my $o = construct_row($type, $self, $row);

        push @{$self->{'track'}}, $o;

        $self->update_fieldwidths($o);
    }
    #vmstat("Panel::append done");
}

sub render_panel {
    my ($self, $par, $posnwidths, $labelwidths) = @_;

    return  unless @{$self->{'track'}};

    #no fieldwidths from caller? use own set
    if (! defined $posnwidths) {
        $posnwidths  = $self->{'posnwidths'};
        $labelwidths = $self->{'labelwidths'};
    }

    if (@{$self->{'header'}}) {
        $par->{'dev'}->render_tr_pre_begin;
        foreach my $s (@{$self->{'header'}}) {
            $par->{'dev'}->render_text($s);
        }
        $par->{'dev'}->render_tr_pre_end;
    }

    $par->{'dev'}->render_tr_pre_begin;
    $self->render_pane($par, $posnwidths, $labelwidths);
    $par->{'dev'}->render_tr_pre_end;

    $self->free_rows;  #garbage collect rows
}

######################################################################
# private methods
######################################################################
sub free_rows {
    #print "free: $_[0]\n";
    while (@{$_[0]->{'track'}}) {  #consume tracks
        my $o = shift @{$_[0]->{'track'}};
        $o = undef;
    }
}

sub initialise_fieldwidths {
    my $self = shift;

    sub swap { return ($_[1], $_[0]) }

    my $par = $self->{'par'};

    #initialise left/right position widths seen so far in this panel,
    #accounting for chunking which creates new positions
    my $chunk      = $self->get_chunk_size;
    my $start      = $self->{'start'};
    my $stop       = $self->{'stop'};
    my $posnwidths = [$MINPOSWIDTH, $MINPOSWIDTH];

    for (my $j=$start; $j < $stop; $j += $chunk) {

        my ($pos0, $pos1) = (length("$j"), length("@{[$j+$chunk-1]}"));

        $posnwidths->[0] = max($posnwidths->[0], $pos0);
        $posnwidths->[1] = max($posnwidths->[1], $pos1);
    }
    @$posnwidths = swap(@$posnwidths)  unless $self->{'forwards'};

    $self->{'chunk'}      = $chunk;
    $self->{'posnwidths'} = $posnwidths;

    #warn "pw: [@{$self->{'posnwidths'}}] fw: $self->{'forwards'}\n";

    #initialise labelwidths
    $self->{'labelwidths'} = [];
    for (my $i=0; $i < @{$par->{'labelflags'}}; $i++) {
        $self->{'labelwidths'}->[$i] = 0;
    }
}

sub update_fieldwidths {
    my ($self, $o) = @_;
    #update label widths seen so far in this panel
    for (my $i=0; $i < @{$self->{'labelwidths'}}; $i++) {
        $self->{'labelwidths'}->[$i] =
            max($self->{'labelwidths'}->[$i], $o->labelwidth($i));
    }
}

sub get_chunk_size {
    my $self = shift;
    my $chunk = $self->{'par'}->{'width'};
    return $self->{'length'}  if $chunk < 1;                     #full width
    return $self->{'length'}  if !$self->{'par'}->{'sequences'}; #sequences off
    return $chunk;
}

#output a pane of chunks
sub render_pane {
    my ($self, $par, $posnwidths, $labelwidths) = @_;

    #need space for sequence position numbers?
    my $has_ruler = 0;
    foreach my $o (@{$self->{'track'}}) {
        $has_ruler = 1  if $o->is_ruler;
        $o->reset;
    }

    #vmstat("render pane");
    while (1) {
        last  unless $self->render_chunk($par, $has_ruler, $self->{'chunk'},
                                         $posnwidths, $labelwidths);
    }
    #vmstat("render pane done");
}

#output a single chunk
sub render_chunk {
    my ($self, $par, $has_ruler, $chunk, $posnwidths, $labelwidths) = @_;

    #render each track's segment for this chunk
    foreach my $o (@{$self->{'track'}}) {

        my $seg = $o->next_segment($par, $chunk);

        return 0  unless defined $seg;  #all chunks done

        #label0: rownum
        if ($par->{'labelflags'}->[0] and $labelwidths->[0]) {
            $par->{'dev'}->render_rownum($labelwidths->[0], $o->label(0));
            $par->{'dev'}->render_hspace($HSPACE);
        }

        #label1: identifier
        if ($par->{'labelflags'}->[1] and $labelwidths->[1]) {
            $par->{'dev'}->render_identifier($labelwidths->[1], $o->label(1),
                                             $o->{'url'});
            $par->{'dev'}->render_hspace($HSPACE);
        }

        #label2: description
        if ($par->{'labelflags'}->[2] and $labelwidths->[2]) {
            $par->{'dev'}->render_description($labelwidths->[2], $o->label(2));
            $par->{'dev'}->render_hspace($HSPACE);
        }

        #labels3-7: info
        for (my $i=3; $i < 8; $i++) {
            if ($par->{'labelflags'}->[$i] and $labelwidths->[$i]) {
                $par->{'dev'}->render_annotation($labelwidths->[$i],
                                                 $o->label($i));
                $par->{'dev'}->render_hspace($HSPACE);
            }
        }

        if ($par->{'sequences'}) {
            #left position
            if ($has_ruler) {
                $par->{'dev'}->render_position($RJUST, $posnwidths->[0],
                                               $seg->[0], $o->is_ruler,
                                               $par->{'bold'})
            }
            $par->{'dev'}->render_hspace($HSPACE);

            #sequence string
            $par->{'dev'}->render_sequence($seg->[2], $par->{'bold'});
            $par->{'dev'}->render_hspace($HSPACE);

            #right position
            if ($has_ruler) {
                $par->{'dev'}->render_position($LJUST, $posnwidths->[1],
                                               $seg->[1], $o->is_ruler,
                                               $par->{'bold'});
            }
        }

        #labels8+: optional info
        for (my $i=8; $i < @$labelwidths; $i++) {
            if ($par->{'labelflags'}->[$i] and $labelwidths->[$i]) {
                $par->{'dev'}->render_hspace($HSPACE);
                $par->{'dev'}->render_text_left($labelwidths->[$i],
                                                $o->label($i));
            }
        }

        $par->{'dev'}->render_newline;
    }

    #blank between chunks
    $par->{'dev'}->render_text("\n");

    return 1;  #end of chunk
}

######################################################################
# private class methods
######################################################################
sub construct_row {
    my ($type, $owner, $data) = @_;
    $type = ucfirst $type;
    no strict 'refs';
    my $row = "Bio::MView::Display::$type"->new($owner, $data);
    use strict 'refs';
    return $row;
}

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

sub dump {
    my $self = shift;
    foreach my $k (sort keys %$self) {
        warn sprintf "%15s => %s\n", $k, $self->{$k};
    }
    map { $_->dump } @{$self->{'track'}};
}

###########################################################################
1;
