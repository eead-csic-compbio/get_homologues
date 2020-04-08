# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Align::Sequence;

use Bio::MView::Align::Row;
use Bio::MView::Color::ColorMap;
use Bio::MView::Align::ColorMixin;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Row Bio::MView::Align::ColorMixin);

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    die "${type}::new: missing arguments\n"  if @_ < 2;
    my ($uid, $sob) = @_;

    my $self = new Bio::MView::Align::Row('sequence', $uid);
    bless $self, $type;

    $self->{'from'}   = $sob->lo;  #start number of sequence
    $self->{'string'} = $sob;      #sequence object
    $self->{'cov'}    = 0;         #percent coverage of reference sequence
    $self->{'pid'}    = 0;         #percent identity to reference sequence

    $self->reset_display;          #hash of display parameters

    $Bio::MView::Align::ColorMixin::FIND_WARNINGS = 1;  #reset

    $self;
}

######################################################################
# public methods
######################################################################
#override
sub is_sequence { 1 }

#overrides Bio::MView::Align::Row::string
sub string   { return $_[0]->{'string'}->string }

sub from     { return $_[0]->{'from'} }
sub seqobj   { return $_[0]->{'string'} }
sub sequence { return $_[0]->{'string'}->sequence }
sub seqlen   { return $_[0]->{'string'}->seqlen }

sub set_coverage {
    my ($self, $ref) = @_;
    my $val = $self->compute_coverage_wrt($ref);
    $self->{'cov'} = $val;
    $self->set_display('label4' => sprintf("%.1f%%", $val));
}

sub get_coverage { return $_[0]->{'cov'} }
sub get_coverage_string { return $_[0]->get_label(4) }

sub set_identity {
    my ($self, $ref, $mode) = @_;
    my $val = $self->compute_identity_to($ref, $mode);
    $self->{'pid'} = $val;
    $self->set_display('label5' => sprintf("%.1f%%", $val));
}

sub get_identity { return $_[0]->{'pid'} }
sub get_identity_string { return $_[0]->get_label(5) }

# Compute the percent coverage of a row with respect to a reference row as:
#   number~of~residues~in~row~aligned~with~reference~row /
#   length~of~ungapped~reference~row
sub compute_coverage_wrt {
    #warn "Bio::MView::Align::Sequence::compute_coverage_wrt(@_)\n";
    my ($self, $othr) = @_;

    return 0      unless defined $othr;
    return 100.0  if $self == $othr;  #always 100% coverage of self

    die "${self}::compute_coverage_wrt: length mismatch\n"
        unless $self->length == $othr->length;

    my $hit = $self->{'string'};
    my $ref = $othr->{'string'};
    my $end = $self->length + 1;

    my ($hitlen, $reflen) = (0, 0);

    for (my $i=1; $i<$end; $i++) {

        #reference must be a sequence character
        my $c = $ref->char_at($i);
        next  unless defined $c;
        $reflen++;

        $c = $hit->char_at($i);
        $hitlen++  if defined $c;
    }

    #compute percent coverage
    return 100.0 * $hitlen/$reflen  if $reflen > 0;
    return 0;
}

#Compute percent identity to a reference row; normalisation depends on the
#mode argument:
#- 'reference' divides by the total reference sequence length
#- 'hit'       divides by the total hit sequence length
#- 'aligned'   divides by the aligned hit region length
sub compute_identity_to {
    #warn "Bio::MView::Align::Sequence::compute_identity_to(@_)\n";
    my ($self, $othr, $mode) = (@_, 'aligned');

    return 0      unless defined $othr;
    return 100.0  if $self == $othr;  #always 100% identical to self

    die "${self}::compute_identity_to: length mismatch\n"
        unless $self->length == $othr->length;

    my $hit = $self->{'string'};
    my $ref = $othr->{'string'};
    my $end = $self->length + 1;

    if ($mode eq 'aligned') {
        return $self->compute_identity_to_region($hit, $ref, $end);
    }
    if ($mode eq 'reference') {
        #note order of arguments
        return $self->compute_identity_to_seq($hit, $ref, $end);
    }
    if ($mode eq 'hit') {
        #note reversed order of arguments
        return $self->compute_identity_to_seq($ref, $hit, $end);
    }
    die "${self}::compute_identity_to: unknown mode '$mode'\n";
}

######################################################################
# private methods
######################################################################
#override: label1 is valid for computed rows
sub reset_display {
    $_[0]->SUPER::reset_display(
        'sequence' => $_[0]->seqobj,
        'label1'   => $_[0]->uid,
    );
}

#override
sub length { return $_[0]->{'string'}->length }

sub compute_identity_to_region {
    my ($self, $hit, $ref, $end) = @_;

    my $ids = 0;  #identities count
    my $len = 0;  #region length
    my $test = $hit;

    for (my $i=1; $i<$end; $i++) {

        my $hitsym = $hit->raw($i);
        my $refsym = $ref->raw($i);

        #ignore pure gaps: at least one sequence character
        my $chars = 0;
        $chars++  if $test->is_char($hitsym);
        $chars++  if $test->is_char($refsym);
        next  unless $chars > 0;

        #length of HIT sequence in aligned region: excludes terminal gaps
        $len++  unless $test->is_terminal_gap($hitsym);

        $hitsym = uc $hitsym; $refsym = uc $refsym;  #standardize case

        #unknown character: contributes to length only
        next  if $hitsym eq 'X' or $refsym eq 'X';

        #count identities
        $ids++  if $hitsym eq $refsym;
    }

    return 100.0 * $ids / $len  if $len > 0;
    return 0;
}

sub compute_identity_to_seq {
    my ($self, $hit, $ref, $end) = @_;

    my $ids = 0;  #identities count
    my $test = $hit;

    for (my $i=1; $i<$end; $i++) {

        my $hitsym = $hit->raw($i);
        my $refsym = $ref->raw($i);

        #ignore pure gaps: at least one sequence character
        my $chars = 0;
        $chars++  if $test->is_char($hitsym);
        $chars++  if $test->is_char($refsym);
        next  unless $chars > 0;

        $hitsym = uc $hitsym; $refsym = uc $refsym;  #standardize case

        #unknown character: contributes to length only
        next  if $hitsym eq 'X' or $refsym eq 'X';

        #count identities
        $ids++  if $hitsym eq $refsym;
    }

    my $len = $ref->seqlen;  #supplied 'ref' sequence

    return 100.0 * $ids / $len  if $len > 0;
    return 0;
}

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

###########################################################################
1;
