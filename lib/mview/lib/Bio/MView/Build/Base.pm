# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::MView::Build::Base;

use Bio::Util::System qw(vmstat);
use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Scheduler;
use Bio::MView::Align::Alignment;

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    die "${type}::new() missing argument\n"  if @_ < 1;
    my $entry = shift;

    my $self = {};
    bless $self, $type;

    $self->{'entry'}     = $entry;  #parse tree ref
    $self->{'align'}     = undef;   #current alignment
    $self->{'index2row'} = undef;   #list of aligned rows, from zero
    $self->{'uid2row'}   = undef;   #hash of aligned rows; by Build::Row->uid
    $self->{'ref_row'}   = undef;   #reference row ref
    $self->{'topn'}      = undef;   #actual number of rows to show
    $self->{'aligned'}   = undef;   #treat input as aligned
    $self->{'keep_uid'}  = undef;   #hash 'keeplist' by Row->uid
    $self->{'nops_uid'}  = undef;   #hash 'nopslist' by Row->uid
    $self->{'hide_uid'}  = undef;   #hash merge 'disc/keep/nops/' by Row->uid

    $self->initialise;

    $self;
}

######################################################################
# public methods
######################################################################
#override if children have a query sequence (children of Build::Search)
sub is_search { 0 }

sub get_entry { $_[0]->{'entry'} }

sub get_ref_row { $_[0]->{'ref_row'} }

sub get_ref_uid {
    return defined $_[0]->{'ref_row'} ? $_[0]->{'ref_row'}->uid : undef
}

sub uid2row   { $_[0]->{uid2row}->{$_[1]} }
sub index2row { $_[0]->{index2row}->[$_[1]] }

sub reset {
    my $self = shift;

    $self->{'aligned'} = 0;

    #how many expected rows of alignment to show (1 more if search)
    $self->{'topn'} = $PAR->get('topn');
    $self->{'topn'} += $self->is_search  if $self->{'topn'} > 0;

    $self->reset_child;
}

#return the block of sequences, 0 if empty block, or undef if no more work
sub next_align {
    my $self = shift;

    #drop old data structures: GC *before* next assignment!
    $self->{'align'} = $self->{'index2row'} = undef;

    #extract an array of row objects
    $self->{'index2row'} = $self->parse;
    #vmstat("Build->next(parse) done");

    #finished?
    return undef  unless defined $self->{'index2row'};

    #$self->dump_index2row;

    #this block empty?
    return 0  unless @{$self->{'index2row'}};

    $self->{'align'} = $self->build_align;
    #vmstat("Build->next(build_align) done");

    #maybe more data but this alignment empty? (identity filtered)
    return 0  unless defined $self->{'align'};
    return 0  unless $self->{'align'}->visible_ids > 0;

    return $self->{'align'};
}

#subclass overrides
sub header {
    my ($self, $quiet) = (@_, 0);
    return ''  if $quiet;

    my $showpcid = $PAR->get('label5');
    my $minident = $PAR->get('minident');
    my $maxident = $PAR->get('maxident');
    my $pcidmode = $PAR->get('pcid');
    my $topn     = $PAR->get('topn');

    my $s = '';

    if (defined $self->{'ref_row'}) {
        my $name = $self->{'ref_row'}->num;
        my $cid  = $self->{'ref_row'}->cid;
        $name = "query"  if $name eq '';
        $s .= "Reference sequence ($name): $cid\n";
    }
    if (0 < $minident and $maxident < 100) {
        $s .= "Identity limits: $minident-$maxident%";
        $s .= " normalised by $pcidmode length.\n";
    } elsif (0 < $minident) {
        $s .= "Minimum identity: $minident%";
        $s .= " normalised by $pcidmode length.\n";
    } elsif ($maxident < 100) {
        $s .= "Maximum identity: $maxident%";
        $s .= " normalised by $pcidmode length.\n";
    } elsif ($showpcid) {
        $s .= "Identities normalised by $pcidmode length.\n";
    }
    if ($topn) {
        $s .= "Maximum sequences to show: $topn\n";
    }

    return $s;
}

#subclass overrides
sub subheader {''}

######################################################################
# protected methods
######################################################################
#return 1 is row should be ignored by row rank or identifier
sub skip_row { my $self = shift; return ! $self->use_row(@_) }

#return 1 if topn rows already generated, 0 otherwise; ignore if if filtering
#on identity; it is assumed the query is implicitly accepted anyway by the
#parser
sub topn_done {
    my ($self, $num) = @_;
    return 0  if $PAR->get('maxident') != 100;
    return 1  if $PAR->get('topn') > 0 and $num > $PAR->get('topn');
    return 0;
}

#subclass overrides: if children need to do something during creation
sub initialise { $_[0]->{scheduler} = new Bio::MView::Build::Scheduler }

#subclass overrides: if children need to do something before each iteration
sub reset_child { $_[0]->{scheduler}->filter }

#subclass overrides: must be overridden
sub use_row { die "$_[0] use_row: virtual method called\n" }

#subclass overrides: must be overridden
sub test_if_aligned { die "$_[0] use_row: virtual method called\n" }

#subclass overrides
sub parse { die "$_[0] parse: virtual method called\n" }

#subclass overrides: map an identifier supplied as {0..N|query|M.N} to
#a list of row objects in $self->{'index2row'}
sub map_id {
    my ($self, $id) = @_;
    my @rowref = ();
    #warn "map_id($id)\n";
    foreach my $row (@{$self->{'index2row'}}) {
        push @rowref, $row  if $row->matches_id($id);
    }
    #warn "${self}::map_id (@rowref)\n";
    return @rowref;
}

#subclass overrides
sub build_rows {
    my ($self, $lo, $hi) = @_;
    foreach my $row (@{$self->{'index2row'}}) {
        unless ($self->{'aligned'}) {
            ($lo, $hi) = $self->get_range($row);
        }
        #warn "Build::build_rows range ($lo, $hi)\n";
        $row->assemble($lo, $hi, $PAR->get('gap'));
    }
}

#subclass overrides
sub get_range {
    my ($self, $row) = @_;
    my @range = @{$PAR->get('range')};
    if (@range and @range % 2 < 1) {
        return ($range[0], $range[1])  if $range[0] < $range[1];
        return ($range[1], $range[0]);
    }
    return $row->range;  #default
}

#subclass overrides: child can return specialised alignment class
sub make_alignment { shift; return new Bio::MView::Align::Alignment(@_) }

######################################################################
# private methods
######################################################################
sub is_aligned { return $_[0]->{'aligned'} == 1 }
sub is_hidden  { return exists $_[0]->{'hide_uid'}->{$_[1]} }
sub is_nop     { return exists $_[0]->{'nops_uid'}->{$_[1]} }

sub outfmt_supported {
    my ($self, $fmt) = @_;
    return 1  if $self->is_aligned;
    return 1  if grep {$_ eq $fmt} qw(fasta pearson pir);
    return 0;
}

#build and return an alignment object using stored rows or undef if no
#alignment data
sub build_align {
    my $self = shift;

    my ($lo, $hi) = $self->get_range($self->{'index2row'}->[0]);

    $self->{'aligned'} = $self->test_if_aligned($lo, $hi);

    #warn "KEEPINSERTS: @{[$PAR->get('keepinserts')]}\n";
    #warn "ALIGNED:     $self->{'aligned'}\n";

    my $outfmt = $PAR->get('outfmt');

    unless ($self->outfmt_supported($outfmt)) {
        warn "Sequence lengths differ for output format '$outfmt' - aborting\n";
        return undef;
    }

    $self->build_indices;
    $self->build_rows($lo, $hi);

    my $aln = $self->make_alignment($self, undef);

    $self->build_base_alignment($aln);

    return undef  unless $aln->size > 0;

    $self->build_mview_alignment($aln)  if $outfmt eq 'mview';

    return $aln;
}

sub build_indices {
    my $self = shift;

    $self->{'uid2row'}  = {};
    $self->{'keep_uid'} = {};
    $self->{'nops_uid'} = {};
    $self->{'hide_uid'} = {};

    #index the row objects by unique 'uid' for fast lookup
    foreach my $i (@{$self->{'index2row'}}) {
        $self->{'uid2row'}->{$i->uid} = $i;
    }

    #make all skiplist rows invisible
    $self->set_uid_fields($PAR->get('skiplist'), 'hide_uid');

    #make all keeplist rows visible again, partially overriding previous
    $self->set_uid_fields($PAR->get('keeplist'), 'keep_uid', 'hide_uid');

    #save any reference row and put it on the keeplist as well
    my $refid = $PAR->get('ref_id');
    if (my @refid = $self->map_id($refid)) {
        $self->{'ref_row'} = $refid[0];  #save the reference row
        $self->set_uid_fields([$refid], 'keep_uid');
    }

    #save the nopslist
    $self->set_uid_fields($PAR->get('nopslist'), 'nops_uid');

    #warn "ref:  ",$self->{'ref_row'}->uid, "\n" if defined $self->{'ref_row'};
    #warn "keep: [", join(",", sort keys %{$self->{'keep_uid'}}), "]\n";
    #warn "nops: [", join(",", sort keys %{$self->{'nops_uid'}}), "]\n";
    #warn "hide: [", join(",", sort keys %{$self->{'hide_uid'}}), "]\n";
}

#iterate over idlist, setting self.finsert for each id; if optional fdelete is
#present remove corresponding self.fdelete element
sub set_uid_fields {
    my ($self, $idlist, $finsert, $fdelete) = @_;
    foreach my $id (@$idlist) {
        my @ids = $self->map_id($id);
        foreach my $r (@ids) {
            my $uid = $r->uid;
            $self->{$finsert}->{$uid} = 1;
            if (defined $fdelete and exists $self->{$fdelete}->{$uid}) {
                delete $self->{$fdelete}->{$uid};
            }
        }
    }
}

sub build_base_alignment {
    my ($self, $aln) = @_;

    #make and insert each sequence into the alignment
    $aln->append_sequences(@{$self->{'index2row'}});

    #filter alignment based on %identity to reference
    $aln->prune_identities($self->get_ref_uid,
                           $PAR->get('pcid'),
                           $PAR->get('minident'),
                           $PAR->get('maxident'),
                           $self->{'topn'},
                           $self->{'keep_uid'});

    my $sequences = $PAR->get('sequences');
    my $noinserts = !$PAR->get('keepinserts');

    #compute columnwise data for aligned output
    if ($self->{'aligned'} and $sequences and $noinserts) {

        if (defined $self->{'ref_row'}) {
            $aln->set_coverage($self->{'ref_row'}->uid);
            $aln->set_identity($self->{'ref_row'}->uid, $PAR->get('pcid'));
        }

        #copy computed data into build row objects
        foreach my $row (@{$self->{'index2row'}}) {
            next  if $self->is_hidden($row->uid);

            my $arow = $aln->uid2row($row->uid);
            next  unless defined $arow;

            $row->set_coverage($arow->get_coverage_string);
            $row->set_identity($arow->get_identity_string);
        }
    }

    #foreach my $r ($aln->all_ids) { $aln->uid2row($r)->seqobj->dump }
    #warn "Alignment width: ", $aln->length;
}

sub build_mview_alignment {
    my ($self, $aln) = @_;

    foreach my $row (@{$self->{'index2row'}}) {
        next  if $self->is_hidden($row->uid);

        my $arow = $aln->uid2row($row->uid);
        next  unless defined $arow;

        my @labels = $row->display_column_values;

        $labels[0] = ''  if $self->is_nop($row->uid);

        #warn "labels: [@{[scalar @labels]}] [@{[join(',', @labels)]}]\n";

        $arow->set_display(
            'label0' => $labels[0],
            'label1' => $labels[1],
            'label2' => $labels[2],
            'label3' => $labels[3],
            'label4' => $labels[4],
            'label5' => $labels[5],
            'label6' => $labels[6],
            'label7' => $labels[7],
            'label8' => $labels[8],
            'url'    => $row->url,
            );
    }
}

######################################################################
# debug
######################################################################
#sub DESTROY { print "destroy: $_[0]\n" }

# sub dump_index2row {
#     my $self = shift;
#     for (my $i=0; $i < @{$self->{'index2row'}}; $i++) {
#         warn "$i:  ", $self->index2row($i)->num, " ",
#             $self->index2row($i)->cid, " [",
#             $self->index2row($i)->seq, "]\n";
#     }
# }

###########################################################################
1;
