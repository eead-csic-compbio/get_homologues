# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: HSSP.pm,v 1.31 2015/06/14 17:09:04 npb Exp $

###########################################################################
package Bio::MView::Build::Row::HSSP;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 0,   1,    'chain',         'chain',      '1S',       '' ],
    ]
}


###########################################################################
package Bio::MView::Build::Format::HSSP;

use Bio::MView::Build::Search;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Search);

#the name of the underlying NPB::Parse::Format parser
sub parser { 'HSSP' }

my %Known_Parameters =
    (
     #name        => [ format  default ]
     'chain'      => [ [],     undef   ],
    );

#tell the parent
sub known_parameters { \%Known_Parameters }

#called by the constructor
sub initialise_child {
    my $self = shift;

    #MaxHom/HSSP chain names
    my @names = $self->{'entry'}->parse(qw(ALIGNMENT))->get_chains;

    #free object
    $self->{'entry'}->free(qw(ALIGNMENT));

    $self->{scheduler} = new Bio::MView::Build::Scheduler(\@names);

    $self;
}

#called on each iteration
sub reset_child {
    my $self = shift;
    #warn "reset_child [@{$self->{'chain'}}]\n";
    $self->{scheduler}->filter($self->{'chain'});
    $self;
}

#current chain being processed
sub chain { $_[0]->{scheduler}->item }

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s  if $quiet;
    $s .= "Chain: " . $self->chain . "\n";
    $s;    
}

sub parse {
    my $self = shift;
    my ($rank, $use, $head, $prot, $align, $match, $id, $seq, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    $head  = $self->{'entry'}->parse(qw(HEADER));
    $prot  = $self->{'entry'}->parse(qw(PROTEIN));
    $align = $self->{'entry'}->parse(qw(ALIGNMENT));

    push @hit, new Bio::MView::Build::Row::HSSP
	(
	 '',
	 $head->{'pdbid'},
	 $head->{'header'},
	 $self->chain,
	);
    $hit[$#hit]->add_frag($align->get_query($self->chain));

    #extract cumulative scores and identifiers from the ranking and
    #corresponding sequences from the already parsed alignment.
    foreach $match (@{$prot->{'ranking'}}) {
	
	$rank++;
	
	if ($match->{'id'} =~ /\|/) {
	    #looks like a genequiz generated HSSP file
	    $id = $match->{'id'};
	} elsif ($match->{'accnum'}) {
	    #looks like a uniprot derived HSSP file
	    $id = "uniprot|$match->{'accnum'}|$match->{'id'}";
	} else {
	    #give up
	    $id = $match->{'id'};
	}

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $id);

	#warn "KEEP: ($rank,$id)\n";

	$seq = $align->get_sequence($rank, $self->chain);

	#skip empty alignments (aligned to different chain)
	next  if $seq =~ /^\s+$/;

	$seq =~ tr/ /-/;    #replace spaces with hyphens

	push @hit, new Bio::MView::Build::Row::HSSP($rank,
                                                    $id,
						    $match->{'protein'},
						    $self->chain);
        $hit[$#hit]->add_frag($seq);
    }
    #map { $_->print } @hit;

    #free objects
    $self->{'entry'}->free(qw(HEADER PROTEIN ALIGNMENT));
    
    return \@hit;
}


###########################################################################
1;
