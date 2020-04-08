# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Align::Consensus;

use Bio::MView::Sequence::Base;
use Bio::MView::Align::Sequence;
use Bio::MView::Color::ColorMap;
use Bio::MView::GroupMap;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Sequence);

#hardwire the consensus line symcolor
my $SYMCOLOR = $Bio::MView::Color::ColorMap::Colour_Black;

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    die "${type}::new: missing arguments\n"  if @_ < 6;
    my ($from, $to, $tally, $group, $thresh, $ignore) = @_;

    if ($thresh < 50 or $thresh > 100) {
        die "${type}::new: threshold '$thresh\%' outside valid range [50..100]\n";
    }

    my $string =
        Bio::MView::GroupMap::consensus($tally, $group, $thresh, $ignore);

    #encode the new "sequence"
    my $sob = new Bio::MView::Sequence;
    $sob->set_find_pad('\.'); $sob->set_pad('.');
    $sob->set_find_gap('\.'); $sob->set_gap('.');
    $sob->insert([$string, $from, $to]);

    my $self = new Bio::MView::Align::Sequence("consensus/$thresh\%", $sob);
    bless $self, $type;

    $self->{'group'}     = $group;   #groupmap
    $self->{'threshold'} = $thresh;  #consensus threshold

    $self;
}

######################################################################
# public methods
######################################################################
#override
sub is_sequence { 0 }

#override
sub is_consensus { 1 }

#override
sub color_by_type {
    my ($self, $kw) = @_;

    my ($color, $end, $i, $cg, @tmp) = ($self->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $SYMCOLOR;

    #warn "color_by_type($self) 1=$kw->{'aln_colormap'} 2=$kw->{'con_colormap'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

        $cg = $self->{'string'}->raw($i);

        #warn "[$i]= $cg\n";

        #white space: no color
        next  if $self->{'string'}->is_space($cg);

        #gap: gapcolour
        if ($self->{'string'}->is_non_char($cg)) {
            push @$color, $i, 'color' => $kw->{'gapcolor'};
            next;
        }

        #use symbol color/wildcard colour
        @tmp = $self->get_color_type($cg,
                                     $kw->{'aln_colormap'},
                                     $kw->{'con_colormap'});

        push @$color, $self->color_tag($kw->{'css1'}, $SYMCOLOR, $i, @tmp);
    }
}

#override
sub color_by_identity {
    my ($self, $kw, $othr) = @_;

    my ($color, $end, $i, $cg, @tmp) = ($self->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $SYMCOLOR;

    #warn "color_by_identity($self, $othr) 1=$kw->{'aln_colormap'} 2=$kw->{'con_colormap'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

       $cg = $self->{'string'}->raw($i);

       #white space: no colour
       next  if $self->{'string'}->is_space($cg);

       #gap: gapcolour
       if ($self->{'string'}->is_non_char($cg)) {
           push @$color, $i, 'color' => $kw->{'gapcolor'};
           next;
       }

       #consensus group symbol is singleton: choose colour
       if ($GROUPMAP->is_singleton($self->{'group'}, $cg)) {

           #refer to reference colormap NOT the consensus colormap
           @tmp = $self->get_color($cg, $kw->{'aln_colormap'});

           push @$color, $self->color_tag($kw->{'css1'}, $SYMCOLOR, $i, @tmp);

           next;
       }

       #symbol not in consensus group: use contrast colour
       push @$color, $i, 'color' => $SYMCOLOR;
    }
}

#this is analogous to Bio::MView::Align::Sequence::color_by_identity()
#but the roles of self (consensus) and other (sequence) are reversed.
sub color_by_consensus_sequence {
    my ($self, $kw, $othr) = @_;

    return  unless $othr;
    return  unless $othr->is_sequence;

    die "${self}::color_by_consensus_sequence: length mismatch\n"
        unless $self->length == $othr->length;

    my ($color, $end, $i, $cg, $cs, $c, @tmp) = ($othr->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $kw->{'symcolor'};

    #warn "color_by_consensus_sequence($self, $othr) 1=$kw->{'aln_colormap'} 2=$kw->{'con_colormap'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

        $cg = $self->{'string'}->raw($i); $cs = $othr->{'string'}->raw($i);

        #warn "[$i]= $cg <=> $cs\n";

        #white space: no colour
        next  if $self->{'string'}->is_space($cs);

        #gap: gapcolour
        if ($self->{'string'}->is_non_char($cs)) {
            push @$color, $i, 'color' => $kw->{'gapcolor'};
            next;
        }

        #symbols in consensus group are stored upcased
        $c = uc $cs;

        #symbol in consensus group: choose colour
        if ($GROUPMAP->in_consensus_group($self->{'group'}, $cg, $c)) {

            #colour by sequence symbol
            @tmp = $self->get_color_consensus_sequence($cs, $cg,
                                                       $kw->{'aln_colormap'},
                                                       $kw->{'con_colormap'});

            push @$color,
                $self->color_tag($kw->{'css1'}, $kw->{'symcolor'}, $i, @tmp);

            next;
        }

        #symbol not in consensus group: use contrast colour
        push @$color, $i, 'color' => $kw->{'symcolor'};
    }
}

#this is analogous to Bio::MView::Align::Sequence::color_by_identity()
#but the roles of self (consensus) and other (sequence) are reversed.
sub color_by_consensus_group {
    my ($self, $kw, $othr) = @_;

    return  unless $othr;
    return  unless $othr->is_sequence;

    die "${self}::color_by_consensus_group: length mismatch\n"
        unless $self->length == $othr->length;

    my ($color, $end, $i, $cg, $cs, $c, @tmp) = ($othr->{'display'}->{'range'});

    push @$color, 1, $self->length, 'color' => $kw->{'symcolor'};

    #warn "color_by_consensus_group($self, $othr) 1=$kw->{'aln_colormap'} 2=$kw->{'con_colormap'}\n";

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

        $cg = $self->{'string'}->raw($i); $cs = $othr->{'string'}->raw($i);

        #warn "[$i]= $cg <=> $cs\n";

        #no sequence symbol: whitespace: no colour
        next  if $self->{'string'}->is_space($cs);

        #gap or frameshift: gapcolour
        if ($self->{'string'}->is_non_char($cs)) {
            push @$color, $i, 'color' => $kw->{'gapcolor'};
            next;
        }

        #symbols in consensus group are stored upcased
        $c = uc $cs;

        #symbol in consensus group: choose colour
        if ($GROUPMAP->in_consensus_group($self->{'group'}, $cg, $c)) {

            #colour by consensus group symbol
            #note: both symbols passed; colormaps swapped
            @tmp = $self->get_color_consensus_group($cs, $cg,
                                                    $kw->{'aln_colormap'},
                                                    $kw->{'con_colormap'});

            push @$color,
                $self->color_tag($kw->{'css1'}, $kw->{'symcolor'}, $i, @tmp);

            next;
        }

        #symbol not in consensus group: use contrast colour
        push @$color, $i, 'color' => $kw->{'symcolor'};
    }
}

######################################################################
# private methods
######################################################################
#colours a row of consensus sequence.
# 1. give consensus symbols their own colour.
# 2. the consensus colormap is just one colour name, use that and ignore CSS.
# 3. the consensus may be a residue name: use prevailing residue colour.
# 4. use the prevailing wildcard residue colour.
# 5. give up.
sub get_color_type {
    my ($self, $c, $mapS, $mapG) = @_;
    #warn "get_color_type($self, $c, $mapS, $mapG)\n";

    return $COLORMAP->get_symbol_color($mapG, $c)
        if $COLORMAP->has_symbol_color($mapG, $c);

    return $COLORMAP->get_palette_color($mapG, 'T')  #ignore CSS setting
        if $COLORMAP->has_palette_color($mapG);

    return $COLORMAP->get_symbol_color($mapS, $c)
        if $COLORMAP->has_symbol_color($mapS, $c);

    return $COLORMAP->get_wildcard_color($mapS)
        if $COLORMAP->has_wildcard_color($mapS);

    return ();  #no match
}

#colours a row of 'normal' sequence only where there is a consensus symbol.
# 1. give residues their own colour.
# 2. use the prevailing wildcard residue colour.
# 3. give up.
sub get_color_consensus_sequence {
    my ($self, $cs, $cg, $mapS, $mapG) = @_;
    #warn "get_color_consensus_sequence($self, $cs, $cg, $mapS, $mapG)\n";

    return $COLORMAP->get_symbol_color($mapS, $cs)
        if $COLORMAP->has_symbol_color($mapS, $cs);

    return $COLORMAP->get_wildcard_color($mapS)
        if $COLORMAP->has_wildcard_color($mapS);

    return ();  #no match
}

#colours a row of 'normal' sequence using colour of consensus symbol.
# 1. give residues the colour of the consensus symbol.
# 2. the consensus may be a residue name: use prevailing residue colour.
# 3. use the prevailing wildcard residue colour.
# 4. give up.
sub get_color_consensus_group {
    my ($self, $cs, $cg, $mapS, $mapG) = @_;
    #warn "get_color_consensus_group($self, $cs, $cg, $mapS, $mapG)\n";

    return $COLORMAP->get_symbol_color($mapG, $cg)
        if $COLORMAP->has_symbol_color($mapG, $cg);

    return $COLORMAP->get_symbol_color($mapS, $cg)
        if $COLORMAP->has_symbol_color($mapS, $cg);

    return $COLORMAP->get_wildcard_color($mapS)
        if $COLORMAP->has_wildcard_color($mapS);

    return ();  #no match
}

###########################################################################
1;
