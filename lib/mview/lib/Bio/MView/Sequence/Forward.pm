# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Sequence::Forward;

use Bio::MView::Sequence::Base;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Sequence);

#overrides
sub is_forwards {1}

#overrides: called once this orientation is known
sub insert { my $self = shift; $self->_insert(@_) }

#overrides
sub raw {
    my ($self, $col) = @_;

    my $p = $self->{'lo'} + $col - 1;

    return ''  if $p < $self->{'lo'} or $p > $self->{'hi'};

    return $Bio::MView::Sequence::Mark_Pad  if
        $p < $self->{'reflo'} or $p > $self->{'refhi'};
    return $self->{'seq'}->{$p}  if exists $self->{'seq'}->{$p};
    return $Bio::MView::Sequence::Mark_Gap;
}

#overrides
sub char_at {
    my ($self, $col) = @_;

    my $p = $self->{'lo'} + $col - 1;
    return undef  unless exists $self->{'seq'}->{$p};

    my $c = $self->{'seq'}->{$p};
    return undef  if $c eq $Bio::MView::Sequence::Mark_Fs1 or $c eq $Bio::MView::Sequence::Mark_Fs2;
    return $c;
}

#overrides
sub col {
    my ($self, $col) = @_;

    my $p = $self->{'lo'} + $col - 1;

#    warn("col(+): $col [",
#        exists $self->{'seq'}->{$col} ? $self->{'seq'}->{$col} : '',
#        "]=> $p [",
#        exists $self->{'seq'}->{$p} ? $self->{'seq'}->{$p} : '',
#        "]\n");

    return ''  if $p < $self->{'lo'} or $p > $self->{'hi'};

    return $self->{'text_pad'}  if
        $p < $self->{'reflo'} or $p > $self->{'refhi'};
    return $self->{'text_gap'}  unless exists $self->{'seq'}->{$p};
    return $self->{'text_spc'}  if
        $self->{'seq'}->{$p} eq $Bio::MView::Sequence::Mark_Spc;
    return $self->{'text_fs1'}  if
        $self->{'seq'}->{$p} eq $Bio::MView::Sequence::Mark_Fs1;
    return $self->{'text_fs2'}  if
        $self->{'seq'}->{$p} eq $Bio::MView::Sequence::Mark_Fs2;
    return $self->{'seq'}->{$p};
}

######################################################################
# protected methods
######################################################################
#overrides: input each argument frag as [string, from, to], where from <= to
sub _insert {
    my $self = shift;
    my ($o, $state) = ('+', 0);

    foreach my $frag (@_) {
        die "insert($o) wrong direction ($frag->[1], $frag->[2])\n"
            if $frag->[1] > $frag->[2];

        my ($len, $lo, $hi, $string);

        if (ref $frag->[-1] and ${ $frag->[-1] } ne '')  {
            $string = ${ $frag->[-1] };
            $len = CORE::length $string;
            $lo = $frag->[1];
            $hi = $lo + $len - 1;
        } else {
            $string = ${ $frag->[0] };
            $len = CORE::length $string;
            $lo = $frag->[1];
            $hi = $frag->[2];
        }
        #warn "frg($o)=@$frag\n";
        #warn "frg($o)=$frag->[1], $frag->[2] => $len [$string]\n";

        $self->_update_sequence_limits($lo, $hi);
        $self->_update_sequence_labels($frag->[3], $frag->[4], 'f1', 't1');
        $self->_update_sequence_labels($frag->[5], $frag->[6], 'f2', 't2');
        $self->_populate_sequence_array($o, $state, $len, $lo, $hi, \$string);

        #warn "insert($o): $self->{'lo'} $self->{'hi'}\n";
    }

    #adjust prefix/suffix positions given new lengths
    $self->{'reflo'} = $self->{'lo'} + $self->{'pfxlen'};
    $self->{'refhi'} = $self->{'hi'} - $self->{'sfxlen'};

    #$self->dump;
}

#overrides
sub _substr {
    my ($self, $start, $stop) = @_;

    return ''  if $start < 1 or $stop < 1;   #negative range args
    return ''  unless $self->{'lo'} > 0;     #empty
    return ''  if $stop  < $self->{'lo'};    #missed (too low)
    return ''  if $start > $self->{'hi'};    #missed (too high)

    $start = $self->{'lo'}  if $start < $self->{'lo'};
    $stop  = $self->{'hi'}  if $stop  > $self->{'hi'};
    #warn "$self _substr(+,lo,hi): $start, $stop";
    #warn "_substr(+,beg,end): $self->{'reflo'}, $self->{'refhi'}";
    $stop++;

    my $s = '';

    for (my $i = $start; $i < $stop; $i++) {
        if ($i < $self->{'reflo'} or $i > $self->{'refhi'}) {
            $s .= $Bio::MView::Sequence::Mark_Pad;
            next;
        }
        #warn "_substr(+): $i [$self->{'seq'}->{$i}]";
        if (exists $self->{'seq'}->{$i}) {
            $s .= $self->{'seq'}->{$i};
        } elsif ($i > $self->{'reflo'} and $i < $self->{'refhi'}) {
            #We never write a gap at the edge. Frameshifts upset the
            #numbering so that [begin,end] may exceed the true range.
            #warn "_substr(+): gap: [$i] ($self->{lo},$self->{hi}), ($self->{reflo},$self->{refhi})\n";
            $s .= $Bio::MView::Sequence::Mark_Gap;
        }
    }
    return $s;
}

#overrides
sub _get_sequence_index {
    my ($self, $lo, $hi, $i) = @_;
    return $lo + $i;
}

###########################################################################
1;
