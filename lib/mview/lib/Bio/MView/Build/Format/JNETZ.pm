# Copyright (C) 1998-2015 Nigel P. Brown
# $Id: JNETZ.pm,v 1.17 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Build::Format::JNETZ;

use Bio::MView::Build::Align;
use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying NPB::Parse::Format parser
sub parser { 'JNETZ' }

sub parse {
    my $self = shift;
    my ($rank, $use, $aln, $id, $seq, $row, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    $aln = $self->{'entry'}->parse(qw(ALIGNMENT));
    
    $rank++; $id  = 'res';
    last  if $self->topn_done($rank);
    next  if $self->skip_row($rank, $rank, $id);
    $seq = $aln->get_query;
    $row = new Bio::MView::Build::Simple_Row($rank, $id, '', $seq);
    #no special subtype: use default
    push @hit, $row;

    $rank++; $id  = 'align';
    last  if $self->topn_done($rank);
    next  if $self->skip_row($rank, $rank, $id);
    $seq = $aln->get_align;
    $row = new Bio::MView::Build::Simple_Row($rank, $id, '', $seq);
    $row->set_subtype('jnet.pred');    #override the default
    push @hit, $row;

    $rank++; $id  = 'conf';
    last  if $self->topn_done($rank);
    next  if $self->skip_row($rank, $rank, $id);
    $seq = $aln->get_conf;
    $row = new Bio::MView::Build::Simple_Row($rank, $id, '', $seq);
    $row->set_subtype('jnet.conf');    #override the default
    push @hit, $row;

    $rank++; $id  = 'final';
    last  if $self->topn_done($rank);
    next  if $self->skip_row($rank, $rank, $id);
    $seq = $aln->get_final;
    $row = new Bio::MView::Build::Simple_Row($rank, $id, '', $seq);
    $row->set_subtype('jnet.pred');    #override the default
    push @hit, $row;

    #map { $_->print } @hit;

    #free objects
    $self->{'entry'}->free(qw(ALIGNMENT));

    return \@hit;
}

#return list of parameters (key, val) pairs special to this parse
sub change_parameters {
    my $self = shift;
    my %opt = ();
    $opt{'label0'} = 0;    #don't report rank
    $opt{'label4'} = 0;    #don't report %identity
    $opt{'label5'} = 0;    #don't report sequence positions
    %opt;
}

#use our own Align subclass instead of the generic one
sub change_alignment_type {
    my ($self, $aln) = @_;
    bless $aln, 'Bio::MView::Build::Format::JNETZ::Align';
    $aln->rebless_align_rows;
    $self;
}

#construct a header string describing this alignment
sub header {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    Bio::MView::Display::displaytext($s);
}


###########################################################################
package Bio::MView::Build::Format::JNETZ::Align::Sequence;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align::Sequence);

sub color_row {
    my $self = shift;
    my %par = @_;

    $par{'css1'}     = 0
	unless defined $par{'css1'};
    $par{'symcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'symcolor'};
    $par{'gapcolor'} = $Bio::MView::Align::Row::Colour_Black
	unless defined $par{'gapcolor'};
    $par{'colormap'} = $Bio::MView::Align::Sequence::Default_Colormap
	unless defined $par{'colormap'};

    my ($color, $end, $i, $c, @tmp) = ($self->{'display'}->{'range'});
    
    push @$color, 1, $self->length, 'color' => $par{'symcolor'};

    for ($end=$self->length+1, $i=1; $i<$end; $i++) {

	$c = $self->{'string'}->raw($i);
	
	#warn "[$i]= $c\n";

	#white space: no color
	next    if $self->{'string'}->is_space($c);

	#gap: gapcolour
	if ($self->{'string'}->is_non_char($c)) {
	    push @$color, $i, 'color' => $par{'gapcolor'};
	    next;
	}
	
	#use symbol color/wildcard colour
	@tmp = $self->get_color($c, $par{'colormap'});
	
	if (@tmp) {
	    if ($par{'css1'}) {
		push @$color, $i, 'class' => $tmp[1];
	    } else {
		push @$color, $i, 'color' => $tmp[0];
	    }
	} else {
	    push @$color, $i, 'color' => $par{'symcolor'};
	}
    }
    
    $self->{'display'}->{'paint'}  = 1;
    $self;
}


###########################################################################
package Bio::MView::Build::Format::JNETZ::Align;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Align);

#the 0 here says don't override any colormap of the same name, to
#allow earler loaded user definitions priority - crude, but it'll do.
Bio::MView::Align::load_colormaps(\*DATA, 0);
#print Bio::MView::Align::list_colormaps(0);

#use our own Align subclass instead of the generic one
sub rebless_align_rows {
    my $self = shift;
    my $i;
    for ($i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	bless $self->{'index2row'}->[$i], 
	    'Bio::MView::Build::Format::JNETZ::Align::Sequence';
    }
    $self;
}

#change the header text
sub header {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s .= "Residues colored by:  property\n";
    $s .= "Structure colored by: type\n";
    Bio::MView::Display::displaytext($s);
}

#ignore generic colouring schemes: use our own
sub set_color_scheme {
    my $self = shift;

    $self->set_parameters(@_);
    
    return $self    if $self->{'coloring'} eq 'none';

    $self->color_by_type('colormap'  => $self->{'colormap'},
			 'symcolor'  => $self->{'symcolor'},
			 'gapcolor'  => $self->{'gapcolor'},
			 'css1'      => $self->{'css1'},
			);

    $self;
}

#propagate colour scheme to row objects
sub color_by_type {
    my $self = shift;
    my $i;
    
    for ($i=0; $i<@{$self->{'index2row'}}; $i++) {
	next unless defined $self->{'index2row'}->[$i];
	next if exists $self->{'nopshash'}->{$self->{'index2row'}->[$i]->id};
	next if exists $self->{'hidehash'}->{$self->{'index2row'}->[$i]->id};
	
	if ($self->{'index2row'}->[$i]->{'type'} eq 'sequence') {
	    #sequence row: use default sequence colours but switch off css
	    $self->{'index2row'}->[$i]->color_row(@_, 'css1'=> 0);
	    next;
	}
	if ($self->{'index2row'}->[$i]->{'type'} eq 'jnet.pred') {
	    #structure row: use our colours
	    $self->{'index2row'}->[$i]->color_row(@_,'colormap'=> 'JNET.PRED');
	    next;
	}
	if ($self->{'index2row'}->[$i]->{'type'} eq 'jnet.conf') {
	    #confidence row: use our colours
	    $self->{'index2row'}->[$i]->color_row(@_,'colormap'=> 'JNET.CONF');
	    next;
	}
    }
    $self;
}


###########################################################################
1;

__DATA__

[JNET.PRED]
Hh  =>  bright-red      #helix
Ee  =>  bright-blue     #sheet
Ll  =>  dark-green      #coil

[JNET.CONF]
0   ->  gray0           #bad
1   ->  gray1
2   ->  gray2
3   ->  gray4
4   ->  gray6
5   ->  gray8
6   ->  gray10
7   ->  gray12
8   ->  gray14
9   ->  gray15          #good


###########################################################################
