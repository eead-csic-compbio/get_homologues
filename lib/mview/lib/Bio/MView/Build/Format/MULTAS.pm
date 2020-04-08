# Copyright (C) 1998-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::MView::Build::Format::MULTAS;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Align;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying Bio::Parse::Format parser
sub parser { 'MULTAS' }

my %Known_Parameters =
    (
     #name        => [ format  default ]
     'block'      => [ [],     undef   ],
    );

#called by the constructor
sub initialise {
    my $self = shift;

    #MULTAL/MULTAS ordinal block number: counted 1..N whereas the actual
    #block has its own 'number' field which is reported by subheader().
    my $last = $self->{'entry'}->count(qw(BLOCK));

    #free object
    $self->{'entry'}->free_parsers(qw(BLOCK));

    $self->{scheduler} = new Bio::MView::Build::Scheduler([1..$last]);

    $self->{'parsed'} = undef;  #ref to parsed block

    $self;
}

#called on each iteration
sub reset_child {
    my $self = shift;
    #warn "reset_child [@{[$PAR->get('block')]}]\n";
    $self->{scheduler}->filter($PAR->get('block'));
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

    $self->{'parsed'} = $self->{'entry'}->parse("BLOCK", $self->block);

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

        push @hit, new Bio::MView::Build::Row::MULTAS($rank+1, $id, $desc,
                                                      $seq);
    }
    #map { $_->dump } @hit;

    #free objects
    $self->{'entry'}->free_parsers(qw(BLOCK));

    return \@hit;
}


###########################################################################
package Bio::MView::Build::Row::MULTAS;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Simple_Row);

sub ignore_columns { ['posn1', 'posn2']; }


###########################################################################
1;
