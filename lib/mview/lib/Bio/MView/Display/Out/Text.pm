# Copyright (C) 2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Display::Out::Text;

use Exporter;
use vars qw(@ISA @EXPORT $LJUST $RJUST);

@ISA = qw(Exporter);

push @EXPORT, qw($LJUST $RJUST);

$LJUST = '-';  #left justify text
$RJUST = '';   #right justify text

sub new {
    my $type = shift;
    die "${type}::new: missing arguments\n"  if @_ < 1;
    my $self = {};
    bless $self, $type;

    $self->{'stm'} = shift;  #o/p stream

    $self;
}

######################################################################
# public methods
######################################################################
#subclass overrides
sub process_char {
    my ($self, $c) = @_;
    return $c;
}

#subclass overrides
sub process_segment {
    my ($self, $string) = @_;
    return $string;
}

###########################################################################
#subclass overrides
sub render_rownum {
    my ($self, $w, $s) = @_;
    #justify right if numeric, left if text
    my $just = ($s =~ /^\d+$/ ? $RJUST : $LJUST);
    $self->write($self->format_label($just, $w, $s));
}

#subclass overrides
sub render_identifier {
    my ($self, $w, $s) = @_;
    $self->render_text_left($w, $s);
}

#subclass overrides
sub render_description {
    my ($self, $w, $s) = @_;
    $self->render_text_left($w, $s);
}

#subclass overrides
sub render_annotation {
    my ($self, $w, $s) = @_;
    $self->render_text_right($w, $s);
}

#subclass overrides
sub render_position {
    my ($self, $just, $w, $s, $isruler) = @_;
    $self->write($self->format_label($just, $w, ($isruler ? $s : '')));
}

#subclass overrides
sub render_sequence {
    my ($self, $s) = @_;
    $self->write($s);
}

#subclass overrides
sub render_hspace {
    my ($self, $w) = @_;
    $self->write(' ' x $w);
}

#subclass overrides
sub render_text_left {
    my ($self, $w, $s) = @_;
    $self->write($self->format_label($LJUST, $w, $s));
}

#subclass overrides
sub render_text_right {
    my ($self, $w, $s) = @_;
    $self->write($self->format_label($RJUST, $w, $s));
}

#subclass overrides: any undef will be rendered as newline
sub render_text {
    my $self = shift;
    foreach my $s (@_) {
        $self->write($s), next  if defined $s;
        $self->write("\n");  #replace undef with newline
    }
}

#subclass overrides
sub render_newline {
    my $self = shift;
    $self->write("\n");
}

#subclass overrides
sub render_table_begin {}

#subclass overrides
sub render_table_end {}

#subclass overrides
sub render_tr_pre_begin {}

#subclass overrides
sub render_tr_pre_end {}

###########################################################################
# protected methods
###########################################################################
sub write {
    my $self = shift;
    my $stm = $self->{'stm'};
    print $stm @_;
}

#left or right justify in padded fieldwidth
sub format_label {
    my ($self, $just, $w, $s) = @_;
    return sprintf("%${just}${w}s", $s);
}

###########################################################################
1;
