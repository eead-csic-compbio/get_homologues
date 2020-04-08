# Copyright (C) 2013-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::MView::Build::Row::MAF;

use Bio::MView::Build::Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub schema {[
    # use? rdb?  key              label         format   default
    [ 1,   1,    'start',         'start',      '8N',       '' ],
    [ 2,   2,    'size',          'size',       '8N',       '' ],
    [ 3,   3,    'strand',        'str',        '3N',       '' ],
    [ 4,   4,    'srcsize',       'srcsize',    '10N',      '' ],
    ]
}

sub ignore_columns { ['desc', 'posn1', 'posn2']; }


###########################################################################
package Bio::MView::Build::Format::MAF;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Align;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Align);

#the name of the underlying Bio::Parse::Format parser
sub parser { 'MAF' }

my %Known_Parameters =
    (
     #name        => [ format  default ]
     'block'      => [ [],     undef   ],
    );

#called by the constructor
sub initialise {
    my $self = shift;

    #MAF ordinal block number: counted 1..N whereas the actual
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
    return $s  if $quiet;
    $s .= "Block: " . $self->block;
    $s .= "  score: $self->{'parsed'}->{'score'}"
        if $self->{'parsed'}->{'score'} ne '';
    $s .= "  pass: $self->{'parsed'}->{'pass'}"
        if $self->{'parsed'}->{'pass'} ne '';
    $s .= "\n";
}

sub parse {
    my $self = shift;
    my ($rank, $use, $desc, $seq, @hit) = (0);

    return  unless defined $self->{scheduler}->next;

    $self->{'parsed'} = $self->{'entry'}->parse("BLOCK", $self->block);

    #block doesn't exist?
    return  unless defined $self->{'parsed'};

    foreach my $row (@{$self->{'parsed'}->{'row'}}) {

        $rank++;

        last  if $self->topn_done($rank);
        next  if $self->skip_row($rank, $rank, $row->{'id'});

        #warn "KEEP: ($rank,$row->{'id'})\n";

        push @hit, new Bio::MView::Build::Row::MAF(
            $rank,
            $row->{'id'},
            '',
            {
                'start'   => $row->{'start'},
                'size'    => $row->{'size'},
                'strand'  => $row->{'strand'},
                'srcsize' => $row->{'srcsize'}
            });
        $hit[$#hit]->add_frag($row->{'seq'});
    }
    #map { $_->dump } @hit;

    #free objects
    $self->{'entry'}->free_parsers(qw(BLOCK));

    return \@hit;
}


###########################################################################
1;
