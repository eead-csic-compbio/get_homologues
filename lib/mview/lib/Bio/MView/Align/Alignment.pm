# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::MView::Align::Alignment;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Align::Ruler;
use Bio::MView::Align::Consensus;
use Bio::MView::Align::Conservation;
use Bio::MView::Align::Sequence;
use Bio::MView::Align::Special;

sub new {
    my $type = shift;
    #warn "${type}::new: @_\n";
    die "${type}::new: missing arguments\n"  if @_ < 2;
    my ($build, $parent) = (shift, shift);

    my $self = {};
    bless $self, $type;

    $self->{'build'}     = $build;   #calling Build object
    $self->{'parent'}    = $parent;  #identifier of parent sequence
    $self->{'uid2index'} = undef;    #hash of identifiers giving row numbers
    $self->{'index2row'} = undef;    #list of aligned rows, from zero
    $self->{'length'}    = 0;        #alignment width
    $self->{'tally'}     = {};       #column tallies for consensus

    #warn "Align.aligned: @{[$build->is_aligned]}\n";
    #warn "Align.parent:  @{[defined $parent ? $parent : 'undef']}\n";

    $self->index_reset;
    $self->append_rows(@_);

    $self;
}

######################################################################
# public methods
######################################################################
sub size { return scalar @{$_[0]->{'index2row'}} }

sub length { return $_[0]->{'length'} }

#make a Sequence from the given Build::Row(s) and append to the alignment
sub append_sequences {
    my $self = shift;
    foreach my $brow (@_) {
        my $arow = $self->make_sequence($brow);
        $self->append_rows($arow);
    }
}

sub uid2row {
    return undef  unless exists $_[0]->{'uid2index'}->{$_[1]};
    return $_[0]->{'index2row'}->[$_[0]->{'uid2index'}->{$_[1]}];
}

sub sort_alignment {
    my ($self, $mode) = @_;

    return  if $mode eq "none";

    my $ref_uid = $self->get_reference_uid;

    return  unless defined $ref_uid;

    my $ref_row = $self->uid2row($ref_uid);

    #warn ">>> ref_uid: ", $ref_uid, "\n";
    #warn ">>> ref_row: ", $ref_row, "\n\n";

    # get all rows except the reference row
    my @unsorted = grep { $_->uid ne $ref_uid } @{$self->{'index2row'}};
    my @sorted;

    # sort all non-reference rows
    if ($mode eq "cov") {
        @sorted = sort {
            #descending
            return -1  if $a->get_coverage() > $b->get_coverage();
            return  1  if $a->get_coverage() < $b->get_coverage();
            #ascending
            return  $a->uid cmp $b->uid;
        }
        @unsorted;
    }
    elsif ($mode eq "pid") {
        @sorted = sort {
            #descending
            return -1  if $a->get_identity() > $b->get_identity();
            return  1  if $a->get_identity() < $b->get_identity();
            #ascending
            return  $a->uid cmp $b->uid;
        }
        @unsorted;
    }
    elsif ($mode eq "cov:pid") {
        @sorted = sort {
            #descending
            return -1  if $a->get_coverage() > $b->get_coverage();
            return  1  if $a->get_coverage() < $b->get_coverage();
            return -1  if $a->get_identity() > $b->get_identity();
            return  1  if $a->get_identity() < $b->get_identity();
            #ascending
            return  $a->uid cmp $b->uid;
        }
        @unsorted;
    }
    elsif ($mode eq "pid:cov") {
        @sorted = sort {
            #descending
            return -1  if $a->get_identity() > $b->get_identity();
            return  1  if $a->get_identity() < $b->get_identity();
            return -1  if $a->get_coverage() > $b->get_coverage();
            return  1  if $a->get_coverage() < $b->get_coverage();
            #ascending
            return  $a->uid cmp $b->uid;
        }
        @unsorted;
    }
    else {
        die "${self}::sort: unknown mode '$mode'\n";
    }

    # prepend the reference row
    unshift @sorted, $ref_row;

    $self->index_rebuild(\@sorted);
}

#return list of all rows as internal ids
sub all_ids {
    my @tmp = ();
    foreach my $r (@{$_[0]->{'index2row'}}) {
        next  unless defined $r;
        push @tmp, $r->uid;
    }
    return @tmp;
}

#return list of visible rows as internal ids, for output
sub visible_ids {
    my @tmp = ();
    foreach my $r (@{$_[0]->{'index2row'}}) {
        next  unless defined $r;
        next  if $_[0]->is_hidden($r->uid);
        next  if $_[0]->is_nop($r->uid);
        push @tmp, $r->uid;
    }
    return @tmp;
}

#ignore id's in remaining arglist
sub set_identity {
    my ($self, $ref, $mode) = @_;
    #warn "Bio::MView::Align::Alignment::set_identity(@_)\n";

    $ref = $self->uid2row($ref);
    return  unless defined $ref;

    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $r->set_identity($ref, $mode);
    }
}

#ignore id's in remaining arglist
sub set_coverage {
    my ($self, $ref) = @_;
    #warn "Bio::MView::Align::Alignment::set_coverage(@_)\n";

    $ref = $self->uid2row($ref);
    return  unless defined $ref;

    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $r->set_coverage($ref);
    }
}

sub header {
    my ($self, $quiet) = @_;
    return ''  if $quiet;
    return ''  unless $PAR->get('html');
    return $self->color_header;
}

sub set_color_scheme {
    my ($self, $ref) = @_;

    my $kw = $PAR->as_dict;

    return  unless $kw->{'html'};

    my $mode = $kw->{'aln_coloring'};

    $self->color_special($kw);

    while (1) {
        $self->color_none($kw),
            last  if $mode eq 'none';

        $self->color_by_type($kw),
            last  if $mode eq 'any';

        $self->color_by_identity($kw, $ref),
            last  if $mode eq 'identity';

        $self->color_by_mismatch($kw, $ref),
            last  if $mode eq 'mismatch';

        $self->color_by_consensus_sequence($kw),
            last  if $mode eq 'consensus';

        $self->color_by_consensus_group($kw),
            last  if $mode eq 'group';

        warn "set_color_scheme: unknown mode '$mode'\n";
        last;
    }

    #find overlays anything else
    $self->color_by_find_block($kw)  if $kw->{'find'} ne '';
}

sub set_consensus_color_scheme {
    my ($self, $aln, $ref) = @_;

    my $kw = $PAR->as_dict;

    return  unless $kw->{'html'};

    my $mode = $kw->{'con_coloring'};

    while (1) {
        last  if $mode eq 'none';

        $self->color_by_type($kw),
            last  if $mode eq 'any';

        $self->color_consensus_by_identity($kw, $aln, $ref),
            last  if $mode eq 'identity';

        warn "set_consensus_color_scheme: unknown mode '$mode'\n";
        last;
    }
}

#return tuple of (orientation, length, start, stop) for the alignment
sub init_display {
    my $self = shift;
    return ()  unless defined $self->{'parent'};

    my $seq = $self->{'parent'}->{'string'};

    return ($seq->is_forwards, $seq->length, $seq->lo, $seq->hi);
}

#Append Row data to the input Display object: done one at a time to reduce
#memory usage instead of accumulating a potentially long list before passing
#to Display::append(), and to permit incremental garbage collection of each
#Align::Row object once it has been appended.
#Garbage collection is enabled by default, unless the optional argument
#$gc_flag is false. This is essential when further processing of Row objects
#will occur, eg., consensus calculations.
sub append_display {
    my ($self, $dis, $gc_flag) = (@_, 1);
    #warn "append_display($dis, $gc_flag)\n";

    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) {
        my $r = $self->{'index2row'}->[$i];

        next  unless defined $r;
        next  if $self->is_hidden($r->uid);  #but let nops through

        #append the row data structure to the Display object
        $dis->append($r->get_display);

        #optional garbage collection
        $self->do_gc($i)  if $gc_flag;
    }
}

sub do_gc {
    my ($self, $i) = @_;
    if (defined $i) {  #just one
        $self->{'index2row'}->[$i] = undef;
        return;
    }
    for (my $i=0; $i<@{$self->{'index2row'}}; $i++) { #all
        $self->{'index2row'}->[$i] = undef;
    }
}

sub prune_identities {
    my ($self, $refid, $mode, $min, $max, $topn, $keep) = @_;

    sub swap { return ($_[1], $_[0]) }

    $min = 0    if $min < 0;
    $max = 100  if $max > 100;

    #trivial cases
    return  unless defined $refid;

    my $ref = $self->uid2row($refid);  #the reference row
    return  unless defined $ref;

    return  unless $min > 0 or $max < 100;

    #still work if limits are back to front
    ($min, $max) = swap($min, $max)  if $min > $max;

    my @obj = ();

    foreach my $row (@{$self->{'index2row'}}) {
        next  unless defined $row;

        #enforce limit on number of rows
        last  if $topn > 0 and @obj == $topn;

        #keep this row regardless of percent identity
        if (exists $keep->{$row->uid}) {
            push @obj, $row;
            next;
        }

        my $pcid = $row->compute_identity_to($ref, $mode);

        #percent identity outside cutoff?
        next  if $pcid < $min or $pcid > $max;

        push @obj, $row;
    }

    $self->index_reset;
    $self->append_rows(@obj);
}

#generate a new alignment with a ruler based on this alignment
sub build_ruler {
    my ($self, $refobj) = @_;
    my $obj = new Bio::MView::Align::Ruler($self->length, $refobj);
    return new Bio::MView::Align::Alignment($self->{'build'},
                                            $self->{'parent'}, $obj);
}

#generate a new alignment using an existing one but with a line of
#clustal-style conservation string (*:.)
sub build_conservation_row {
    my $self = shift;

    my $moltype = $PAR->get('moltype');

    #extract sequence rows
    my @ids = $self->visible_computable_ids;

    my $from = $self->{'parent'}->from;
    my $to   = $from + $self->{'length'} - 1;
    #warn "fm/to: ($from, $to)\n";

    #alignment column numbering
    my $string = $self->conservation(\@ids, 1, $self->{'length'}, $moltype);

    #sequence object lo/hi numbering
    my $obj = new Bio::MView::Align::Conservation($from, $to, $string);

    return new Bio::MView::Align::Alignment($self->{'build'},
                                            $self->{'parent'}, $obj);
}

#generate a new alignment using an existing one but with lines showing
#consensus sequences at specified percent thresholds
sub build_consensus_rows {
    my $self = shift;

    my $group     = $PAR->get('con_groupmap');
    my $threshold = $PAR->get('con_threshold');
    my $ignore    = $PAR->get('con_ignore');
    my $con_gaps  = $PAR->get('con_gaps');

    my $tally = $self->compute_tallies($PAR->get('con_groupmap'));

    my @obj = ();

    my $from = $self->{'parent'}->from;
    my $to   = $from + $self->length - 1;

    foreach my $thresh (@$threshold) {
        push @obj,new Bio::MView::Align::Consensus($from, $to, $tally,
                                                   $group, $thresh, $ignore);
    }

    return new Bio::MView::Align::Alignment($self->{'build'},
                                            $self->{'parent'}, @obj);
}

sub conservation {
    #conservation mechanism inspired by Clustalx-2.1/AlignmentOutput.cpp
    my $CONS_STRONG = [ qw(STA NEQK NHQK NDEQ QHRK MILV MILF HY FYW) ];
    my $CONS_WEAK   = [ qw(CSA ATV SAG STNK STPA SGND SNDEQK NDEQHK
                           NEQHRK FVLIM HFY) ];

    #(from,to) must be 1-based along alignment
    my ($self, $ids, $from, $to, $moltype) = @_;

    return ''  unless @$ids; #empty alignment

    my @tmp = $self->visible_computable_ids;

    my $refseq = $self->uid2row($tmp[0])->seqobj;
    my $depth = scalar @tmp;
    my $s = '';
    #warn "conservation: from=$from, to=$to, depth=$depth\n";

    my $initcons = sub {
        my $values = shift;
        my $dict = {};
        for my $group (@$values) { $dict->{$group} = 0 }
        $dict;
    };

    my $addcons = sub {
        my ($dict, $char) = @_;
        for my $group (keys %$dict) {
            $dict->{$group}++  if index($group, $char) > -1;
        }
        $dict;
    };

    my $testcons = sub {
        my ($dict, $max) = @_;
        for my $group (keys %$dict) {
            return 1  if $dict->{$group} == $max;
        }
        return 0;
    };

    my $printcons = sub {
        my ($j, $dict, $name, $stm) = (@_, \*STDOUT);
        print "$j, $name\n";
        for my $group (sort keys %$dict) {
            printf $stm "%-6s => %d\n", $group, $dict->{$group};
        }
        print "\n\n";
    };

    #iterate over alignment columns
    for (my $j=$from; $j<=$to; $j++) {

        last  if $j > $refseq->length;

        my $strong  = &$initcons($CONS_STRONG);
        my $weak    = &$initcons($CONS_WEAK);
        my $refchar = $refseq->raw($j);
        my $same    = 0;

        #iterate over sequence list
        for (my $i=0; $i<@$ids; $i++) {
            my $thischar = uc $self->uid2row($ids->[$i])->seqobj->raw($j);
            #warn "[$j][$i] $refchar, $thischar, $ids->[$i]\n";
            next  if $self->is_nop($ids->[$i]);
            $same++   if $thischar eq $refchar;
            &$addcons($strong, $thischar);
            &$addcons($weak, $thischar);
        }
        #&$printcons($j, $strong, 'strong');
        #&$printcons($j, $weak, 'weak');

        #warn "$same, $depth, [@{[$PAR->get('ref_id')]}], $refchar, ", $refseq->is_char($refchar), "\n";
        if ($depth > 0) {
            if ($same == $depth and $refseq->is_char($refchar)) {
                $s .= '*';
                next;
            }
            if ($moltype eq 'aa') {
                $s .= ':', next  if &$testcons($strong, $depth);
                $s .= '.', next  if &$testcons($weak, $depth);
            }
        }
        $s .= ' ';
    }
    #warn "@{[length($s)]} [$s]\n";
    return \$s;
}

######################################################################
# protected methods
######################################################################
#subclass overrides: sequence factory
sub make_sequence {
    my ($self, $row) = @_;

    my $rid = $row->rid;
    my $uid = $row->uid;
    my $sob = $row->sob;

    my $subtype = '';
    $subtype = 'special'  if $rid =~ /^\#/; #leading hash: special row

    #warn "make_sequence: $rid => @{[$subtype ne '' ? $subtype : \"''\"]}\n";

    if ($subtype eq '') {
        return new Bio::MView::Align::Sequence($uid, $sob);
    }

    if ($subtype ne '') {
        return new Bio::MView::Align::Special($uid, $sob);
    }

    die "${self}::make_sequence: unknown type '$subtype'\n";
}

######################################################################
# private methods
######################################################################
sub is_aligned { return $_[0]->{'build'}->is_aligned }
sub is_hidden  { return $_[0]->{'build'}->is_hidden($_[1]) }
sub is_nop     { return $_[0]->{'build'}->is_nop($_[1]) }

sub index_reset {
    $_[0]->{'uid2index'} = {};
    $_[0]->{'index2row'} = [];
}

sub index_insert_row {
    my ($self, $i, $row) = @_;
    $self->{'uid2index'}->{$row->uid} = $i;
    $self->{'index2row'}->[$i] = $row;
}

sub index_rebuild {
    my ($self, $rows) = @_;

    $self->index_reset;

    for (my $i = 0; $i < @$rows; $i++) {
        my $row = $rows->[$i];
        $self->index_insert_row($i, $row);
    }
}

sub append_rows {
    my $self = shift;
    foreach my $row (@_) {
        my $i = @{$self->{'index2row'}};  #number of rows so far
        #warn "append_rows: [$i]\t@{[$row->uid]}\t@{[$row->string]}\n";

        if ($i < 1) {  #empty
            $self->{'length'} = $row->length;
            $self->{'parent'} = $row  unless defined $self->{'parent'};
        }

        if ($self->is_aligned and $row->length != $self->{'length'}) {
            die "${self}::append_rows: incompatible alignment lengths, row $i, expected $self->{'length'}, got @{[$row->length]}\n";
        }

        $self->index_insert_row($i, $row);
    }
}

sub get_reference_uid {
    my $row = $_[0]->{'build'}->get_ref_row();
    return undef  unless defined $row;
    return $row->uid;
}

#return list of visible rows as row objects, for output
sub visible_rows {
    my @tmp = ();
    foreach my $r (@{$_[0]->{'index2row'}}) {
        next  unless defined $r;
        next  if $_[0]->is_hidden($r->uid);
        next  if $_[0]->is_nop($r->uid);
        push @tmp, $r;
    }
    return @tmp;
}

#return list of visible sequence rows excluding nops as internal ids,
#for processing and output
sub visible_computable_ids {
    my @tmp = ();
    foreach my $r (@{$_[0]->{'index2row'}}) {
        next  unless defined $r;
        next  unless $r->is_sequence;
        next  if $_[0]->is_hidden($r->uid);
        next  if $_[0]->is_nop($r->uid);
        push @tmp, $r->uid;
    }
    return @tmp;
}

sub color_header {
    my $self = shift;

    my $mode      = $PAR->get('aln_coloring');
    my $threshold = $PAR->get('aln_threshold');
    my $ref_id    = $PAR->get('ref_id');
    my $find      = $PAR->get('find');

    my $s = '';

    if ($mode eq 'any') {
        $s .= "Colored by: property";
    }
    elsif ($mode eq 'identity' and defined $ref_id) {
        $s .= "Colored by: identity";
    }
    elsif ($mode eq 'mismatch' and defined $ref_id) {
        $s .= "Colored by: mismatch";
    }
    elsif ($mode eq 'consensus') {
        $s .= "Colored by: consensus/$threshold\%";
    }
    elsif ($mode eq 'group') {
        $s .= "Colored by: consensus group/$threshold\%";
    }

    #append any find pattern colouring
    if ($find ne '') {
        $s .= $s eq '' ? "Colored by: " : "; ";
        $s .= "search pattern '$find'";
    }
    $s .= "\n"  if $s ne '';

    return $s;
}

#propagate colour scheme to row objects
sub color_special {
    my ($self, $kw) = @_;
    foreach my $r ($self->visible_rows) {
        next  unless $r->is_special;
        $r->color_special($kw);
    }
}

#propagate colour scheme to row objects
sub color_none {
    my ($self, $kw) = @_;
    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $r->color_none($kw);
    }
}

#propagate colour scheme to row objects
sub color_by_type {
    my ($self, $kw) = @_;
    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence or $r->is_consensus;
        $r->color_by_type($kw);
    }
}

#propagate colour scheme to row objects
sub color_by_identity {
    my ($self, $kw, $id) = @_;

    my $ref = $self->uid2row($id);
    return  unless defined $ref;

    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $r->color_by_identity($kw, $ref);
    }
}

#propagate colour scheme to row objects
sub color_by_mismatch {
    my ($self, $kw, $id) = @_;

    my $ref = $self->uid2row($id);
    return  unless defined $ref;

    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $r->color_by_mismatch($kw, $ref);
    }
}

#propagate colour scheme to row objects
sub color_by_consensus_sequence {
    my ($self, $kw) = @_;

    my $tally = $self->compute_tallies($kw->{'aln_groupmap'});

    my $from = $self->{'parent'}->from;
    my $to   = $from + $self->length - 1;

    my $con = new Bio::MView::Align::Consensus($from, $to, $tally,
                                               $kw->{'aln_groupmap'},
                                               $kw->{'aln_threshold'},
                                               $kw->{'aln_ignore'});
    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $con->color_by_consensus_sequence($kw, $r);
    }
}

#propagate colour scheme to row objects
sub color_by_consensus_group {
    my ($self, $kw) = @_;

    my $tally = $self->compute_tallies($kw->{'aln_groupmap'});

    my $from = $self->{'parent'}->from;
    my $to   = $from + $self->length - 1;

    my $con = new Bio::MView::Align::Consensus($from, $to, $tally,
                                               $kw->{'aln_groupmap'},
                                               $kw->{'aln_threshold'},
                                               $kw->{'aln_ignore'});
    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $con->color_by_consensus_group($kw, $r);
    }
}

#propagate colour scheme to row objects
sub color_by_find_block {
    my ($self, $kw) = @_;
    foreach my $r ($self->visible_rows) {
        next  unless $r->is_sequence;
        $r->color_by_find_block($kw);
    }
}

#propagate colour scheme to row objects
sub color_consensus_by_identity {
    my ($self, $kw, $aln, $id) = @_;

    my $ref = $aln->uid2row($id);
    return  unless defined $ref;

    foreach my $r ($self->visible_rows) {
        next  unless $r->is_consensus;
        $r->color_by_identity($kw, $ref);
    }
}

sub compute_tallies {
    my ($self, $gname) = @_;

    my $gaps = $PAR->get('con_gaps');

    my $key = join('::', ($gname, $gaps));
    #warn "Tally: $key\n";

    #is there already a suitable tally?
    return $self->{'tally'}->{$key}  if exists $self->{'tally'}->{$key};

    $self->{'tally'}->{$key} = [];

    #iterate over columns
    for (my $c=1; $c <= $self->{'length'}; $c++) {

        my $column = [];

        #iterate over rows
        foreach my $r ($self->visible_rows) {
            next  unless $r->is_sequence;
            push @$column, $r->{'string'}->raw($c);
        }

        #warn "compute_tallies: @$column\n";

        push @{$self->{'tally'}->{$key}},
            Bio::MView::GroupMap::tally_column($gname, $column, $gaps);
    }

    return $self->{'tally'}->{$key};
}

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

sub dump {
    sub _format {
        my ($self, $k, $v) = @_;
        $v = 'undef' unless defined $v;
        $v = "'$v'" if $v =~ /^\s*$/;
        return sprintf("  %-15s => %s\n", $k, $v)
    }
    my $self = shift;
    warn "$self\n";
    map { warn $self->_format($_, $self->{$_}) } sort keys %{$self};
    warn "  sequence rows:\n";
    for (my $i=0; $i < @{$self->{'index2row'}}; $i++) {
        my $r = $self->{'index2row'}->[$i];
        print STDERR "  $i: ";
        if (defined $r) {
            print $r->dump;
        } else {
            print STDERR "\n";
        }
    }
    warn "\n";
}

######################################################################
1;
