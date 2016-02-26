# Copyright (C) 1998-2015 Nigel P. Brown
# $Id: MULTAS.pm,v 1.12 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Build::Format::MULTAS;

use Bio::MView::Build::Align;
use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying NPB::Parse::Format parser
sub parser { 'MULTAS' }

my %Known_Parameters =
    (
     #name        => [ format  default ]
     'block'      => [ [],     undef   ],
    );

#tell the parent
sub known_parameters { \%Known_Parameters }

#called by the constructor
sub initialise_child {
    my $self = shift;

    #MULTAL/MULTAS ordinal block number: counted 1..N whereas the actual
    #block has its own 'number' field which is reported by subheader().
    my $last = $self->{'entry'}->count(qw(BLOCK));

    #free object
    $self->{'entry'}->free(qw(BLOCK));

    $self->{scheduler} = new Bio::MView::Build::Scheduler([1..$last]);

    $self->{'parsed'} = undef;  #ref to parsed block

    $self;
}

#called on each iteration
sub reset_child {
    my $self = shift;
    #warn "reset_child [@{$self->{'block'}}]\n";
    $self->{scheduler}->filter($self->{'block'});
    $self;
}

#current block being processed
sub block { $_[0]->{scheduler}->item }

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s .= "Block: " . $self->block . " (rank=$self->{'parsed'}->{'number'})\n";
    $s;    
}

sub parse {
    my $self = shift;
    my ($rank, $use, $list, $align, $id, $desc, $seq, @hit) = ();

    return  unless defined $self->{scheduler}->next;

    $self->{'parsed'} = $self->{'entry'}->parse("BLOCK[@{[$self->block]}]");

    #block doesn't exist?
    return  unless defined $self->{'parsed'};

    $list  = $self->{'parsed'}->parse(qw(LIST));
    $align = $self->{'parsed'}->parse(qw(ALIGNMENT));
    
    if ($list->{'count'} != $align->{'count'}) {
	die "${self}::parser() different alignment and identifier counts\n";
    }
    
    for ($rank=0; $rank < $list->{'count'}; $rank++) {
	
	$id   = $list->{'hit'}->[$rank]->{'id'};
	$desc = $list->{'hit'}->[$rank]->{'desc'};
	$seq  = $align->{'seq'}->[$rank];

        last  if $self->topn_done($rank+1);
        next  if $self->skip_row($rank+1, $rank+1, $id);

	#warn "KEEP: ($rank,$id)\n";

	push @hit, new Bio::MView::Build::Simple_Row($rank+1, $id, $desc, $seq);
    }
    #map { $_->print } @hit;

    #free objects
    $self->{'entry'}->free(qw(BLOCK));

    return \@hit;
}


###########################################################################
1;
