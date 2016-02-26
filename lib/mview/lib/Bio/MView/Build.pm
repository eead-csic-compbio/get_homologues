# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Build.pm,v 1.24 2015/06/18 21:26:11 npb Exp $

######################################################################
package Bio::MView::Build;

use Universal;
use NPB::Parse::Regexps;
use NPB::Parse::Stream;
use Bio::MView::Align;
use Bio::MView::Display;
use Bio::MView::Build::Scheduler;
use strict;

my $DEBUG_PARAM = 0;

my %Template = 
    (
     'entry'       => undef,   #parse tree ref
     'align'       => undef,   #current alignment
     'index2row'   => undef,   #list of aligned rows, from zero
     'uid2row'     => undef,   #hash of aligned rows, by Build::Row->uid
     'ref_row'     => undef,   #reference row ref
     'topn'        => undef,   #show at most top N items
     'show'        => undef,   #actual number of rows to show
     'minident'    => undef,   #show items with at least minident %identity 
     'maxident'    => undef,   #show items with at most maxident %identity 
     'pcid'        => undef,   #identity calculation method
     'mode'        => undef,   #display format mode
     'aligned'     => undef,   #treat input as aligned
     'ref_id'      => undef,   #reference id for %identity

     'skiplist'    => undef,   #discard rows by {num,id,regex}
     'keeplist'    => undef,   #keep rows by {num,id,regex}
     'nopslist'    => undef,   #no-process rows by {num,id,regex}

     'keep_uid'    => undef,   #hashed version of 'keeplist' by Row->uid
     'nops_uid'    => undef,   #hashed version of 'nopslist'  by Row->uid
     'hide_uid'    => undef,   #hashed merge of 'disc/keep/nops/' by Row->uid

     'range'       => undef,   #display lower/upper bounds (sequence numbering)
     'gap'         => undef,   #output sequence gap character
    );

my %Known_Parameters =
    (
     #name        => [ format,     default   ]
     'topn'       => [ '\d+',      0         ],
     'minident'   => [ $RX_Ureal,  0         ],
     'maxident'   => [ $RX_Ureal,  100       ],
     'pcid'       => [ '\S+',      'aligned' ],
     'mode'       => [ '\S+',      'new'     ],
     'ref_id'     => [ '\S+',      0         ],
     'skiplist'   => [ [],         []        ],
     'keeplist'   => [ [],         []        ],
     'nopslist'   => [ [],         []        ],
     'range'      => [ [],         []        ],
     'gap'        => [ '\S',       '+'       ],
     'showpcid'   => [ '\d+',      1,        ],
    );

my %Known_Identity_Mode =
    (
     'reference'  => 1,
     'aligned'    => 1,
     'hit'        => 1,
    );

my %Known_HSP_Tiling =
    (
     'all'       => 1,
     'ranked'    => 1,
     'discrete'  => 1,
    );

sub new {
    my $type = shift;
    #warn "${type}::new(@_)\n";
    if (@_ < 1) {
	die "${type}::new() missing argument\n";
    }
    my $self = { %Template };

    $self->{'entry'} = shift;

    bless $self, $type;

    $self->initialise_parameters;
    $self->initialise_child;

    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

sub print {
    my $self = shift;
    local $_;
    foreach (sort keys %$self) {
	printf "%15s => %s\n", $_, $self->{$_};
    }
    print "\n";
    $self;
}

sub check_identity_mode {
    if (defined $_[0]) {
	if (exists $Known_Identity_Mode{$_[0]}) {
	    return lc $_[0];
	}
    }
    return map { lc $_ } sort keys %Known_Identity_Mode;
}

sub check_hsp_tiling {
    if (defined $_[0]) {
	if (exists $Known_HSP_Tiling{$_[0]}) {
	    return lc $_[0];
	}
    }
    return map { lc $_ } sort keys %Known_HSP_Tiling;
}

sub initialise_parameters {
    my $self = shift;
    foreach my $p (\%Known_Parameters, $self->known_parameters) {
        warn "initialise_parameters $p\n"  if $DEBUG_PARAM;
        $self->_initialise_parameters($p);
    }
    $self;
}

sub _initialise_parameters {
    my ($self, $p) = @_;
    foreach my $k (keys %$p) {
	warn "_initialise_parameters() $k\n"  if $DEBUG_PARAM;
	if (ref $p->{$k}->[0] eq 'ARRAY') {
	    $self->{$k} = [];
	    next;
	}
	if (ref $p->{$k}->[0] eq 'HASH') {
	    $self->{$k} = {};
	    next;
	}
	$self->{$k} = $p->{$k}->[1];
    }
    $self;
}

sub set_parameters {
    my $self = shift;
    foreach my $p (\%Known_Parameters, $self->known_parameters) {
        warn "set_parameters: $p\n"  if $DEBUG_PARAM;
        $self->_set_parameters($p, @_);
    }

    $self->{'aligned'} = 0;

    #how many expected rows of alignment
    $self->{'show'} = $self->{'topn'} + $self->has_query;
    $self;
}

sub _set_parameters {
    my ($self, $p) = (shift, shift);
    while (my $key = shift) {
	my $val = shift;
	if (exists $p->{$key}) {
	    warn "_set_parameters() $key, @{[defined $val ? $val : 'undef']}\n"
                if $DEBUG_PARAM;
	    if (ref $p->{$key}->[0] eq 'ARRAY' and ref $val eq 'ARRAY') {
		$self->{$key} = $val;
		next;
	    }
	    if (ref $p->{$key}->[0] eq 'HASH' and ref $val eq 'HASH') {
		$self->{$key} = $val;
		next;
	    }
	    if (! defined $val) {
		#set default
		$self->{$key} = $p->{$key}->[1];
		next;
	    }
	    if ($val =~ /^$p->{$key}->[0]$/) {
		#matches expected format
		$self->{$key} = $val;
		next;
	    }
	    warn "${self}::set_parameters() bad value for '$key', got '$val', wanted '$p->{$key}->[0]'\n";
	}
	#ignore unrecognised parameters which may be recognised by subclass
	#set_parameters() methods.
    }
    $self;
}

#override if children add extra parameters
sub known_parameters {{}}

#override if children have an extra query sequence (children of Build::Search)
sub has_query {0}

#override if children need to do something special after during creation
sub initialise_child {
    $_[0]->{scheduler} = new Bio::MView::Build::Scheduler;
    $_[0];
}

#override if children need to do something special before each iteration
sub reset_child {
    $_[0]->{scheduler}->filter;
    $_[0];
}

#return 1 if topn rows already generated, 0 otherwise; ignore if if filtering
#on identity; it is assumed the query is implicitly accepted anyway by the
#parser
sub topn_done {
    my ($self, $num) = @_;
    return 0  if $self->{'maxident'} != 100;
    return 1  if $self->{'topn'} > 0 and $num > $self->{'topn'};
    return 0;
}

#return 1 is row should be ignored by row rank or identifier
sub skip_row { my $self = shift; ! $self->use_row(@_) }

#override in children
sub use_row { die "$_[0] use_row() virtual method called\n" }

#map an identifier supplied as {0..N|query|M.N} to a list of row objects in
#$self->{'index2row'}
sub map_id {
    my ($self, $ref) = @_;
    my ($i, @rowref) = ();

    #warn "map_id($ref)\n";

    for ($i=0; $i<@{$self->{'index2row'}}; $i++) {
	
	#major row number = query
	if ($ref =~ /^0$/) {
	    if ($self->{'index2row'}->[$i]->num eq '' or
		$self->{'index2row'}->[$i]->num eq $ref) {
		push @rowref, $self->{'index2row'}->[$i];
	    }
	    next;
	}
	
	#major row number
	if ($ref =~ /^\d+$/) {
	    #exact match
	    if ($self->{'index2row'}->[$i]->num eq $ref) {
		push @rowref, $self->{'index2row'}->[$i];
		next;
	    }
	    #match to major.minor prefix
	    if ($self->{'index2row'}->[$i]->num =~ /^$ref\./) {
		push @rowref, $self->{'index2row'}->[$i];
		next;
	    }
	    next;
	}
	
	#major.minor row number
	if ($ref =~ /^\d+\.\d+$/) {
	    if ($self->{'index2row'}->[$i]->num eq $ref) {
		push @rowref, $self->{'index2row'}->[$i];
	    }
	    next;
	}
	
	#string identifier
	if ($ref eq $self->{'index2row'}->[$i]->rid or 
	    $ref eq $self->{'index2row'}->[$i]->cid) {
	    push @rowref, $self->{'index2row'}->[$i];
	    next;
	}
	
	#regex inside // pair, applied case-insensitive
	if ($ref =~ /^\/.*\/$/) {
	    my $r = $ref;
	    $r =~ s/^\///; $r =~ s/\/$//;
	    if ($self->{'index2row'}->[$i]->cid =~ /$r/i) {
		#warn "map_id: [$i] /$r/ @{[$self->{'index2row'}->[$i]->cid]}\n";
		push @rowref, $self->{'index2row'}->[$i];
	    }
	    next;
	}

	#wildcard
	if ($ref =~ /^\*$/ or $ref =~ /^all$/i) {
	    push @rowref, $self->{'index2row'}->[$i];
	    next;
	}
	
    }
    #warn "${self}::map_id (@rowref)\n";
    return @rowref;
}

#override this to return list of parameters (key, val) pairs special to
#a particular instance (eg., specific to some parser)
sub change_parameters {()}

#allow instance to rebless an Bio::MView::Align object
sub change_alignment_type {}

sub get_entry { $_[0]->{'entry'} }

sub get_row_id {
    my ($self, $id) = @_;
    if (defined $id) {
	my @id = $self->map_id($id);
	return undef        unless @id;
	return $id[0]->uid  unless wantarray;
	return map { $_->uid } @id;
    }
    return undef;
}

sub get_row {
    my ($self, $id) = @_;
    if (defined $id) {
	my @id = $self->map_id($id);
	return undef   unless @id;
	return $id[0]  unless wantarray;
	return @id;
    }
    return undef;
}

sub uid2row   { $_[0]->{uid2row}->{$_[1]} }
sub index2row { $_[0]->{index2row}->[$_[1]] }

#construct a header string describing this alignment
sub header {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    if (defined $self->{'ref_row'}) {
	$s .= "Reference sequence ";
	if ($self->{'ref_row'}->num !~ /^\s*$/) {
	    $s .= "(" . $self->{'ref_row'}->num . ")";
	} else {
	    $s .= "(query)";
	}
	$s .= ": " . $self->{'ref_row'}->cid . "\n";
    }
    if ($self->{'minident'} > 0 and $self->{'maxident'} < 100) {
	$s .= "Pairwise identity limits: $self->{'minident'}-$self->{'maxident'}%";
	$s .= " normalised by $self->{'pcid'} length.\n";
    } elsif ($self->{'minident'} > 0) {
	$s .= "Minimum pairwise identity: $self->{'minident'}%";
	$s .= " normalised by $self->{'pcid'} length.\n";
    } elsif ($self->{'maxident'} < 100) {
	$s .= "Maximum pairwise identity: $self->{'maxident'}%";
	$s .= " normalised by $self->{'pcid'} length.\n";
    } elsif ($self->{'showpcid'}) {
	$s .= "Identities normalised by $self->{'pcid'} length.\n";
    }
    if ($self->{'topn'}) {
	$s .= "Maximum sequences to show: $self->{'topn'}\n";
    }
    Bio::MView::Display::displaytext($s);
}

sub initialise {
    my ($self, $params) = @_;
    $self->set_parameters(%$params);
    $self->reset_child;
}

sub subheader {''}

#return the block of sequences, 0 if empty block, or undef if no more work
sub next {
    my $self = shift;

    #drop old data structures: GC *before* next assignment!
    $self->{'align'} = $self->{'index2row'} = undef;
    
    #extract an array of row objects
    $self->{'index2row'} = $self->parse;
    #Universal::vmstat("Build->next(parse) done");

    #finished? note: "$self->{'align'}->free" is not needed
    return undef  unless defined $self->{'index2row'};

    #my $i; for ($i=0; $i < @{$self->{'index2row'}}; $i++) {
    #    warn "[$i]  ", $self->index2row($i)->num, " ",
    #	      $self->index2row($i)->cid, "\n";
    #}

    #this block empty?
    return 0  unless @{$self->{'index2row'}};

    $self->{'align'} = $self->build_block;
    #Universal::vmstat("Build->next(build_block) done");

    #maybe more data but this alignment empty? (identity filtered)
    return 0  unless defined $self->{'align'};
    return 0  unless $self->{'align'}->visible_ids > 0;

    return $self->{'align'};
}

sub build_block {
    my $self = shift;

    my $aligned = 1;

    my ($lo, $hi) = $self->get_range($self->{'index2row'}->[0]);

    #not a search, so do all rows have same range?
    if ($self->isa('Bio::MView::Build::Align')) {
        for (my $i=1; $i < @{$self->{'index2row'}}; $i++) {
            my ($lo2, $hi2) = $self->get_range($self->{'index2row'}->[$i]);
            #warn "$self->{'index2row'}->[$i] ($lo2, $hi2)\n";
            $aligned = 0, last  if $lo != $lo2 or $hi != $hi2;
        }
    }

    $self->{'aligned'} = $aligned;

    if (!$aligned and
        grep {$_ eq $self->{'mode'}} qw(new clustal msf rdb plain)) {
        warn "Output format is '$self->{'mode'}', but sequence lengths differ - aborting\n";
        return undef;
    }

    $self->build_indices;
    $self->build_rows($lo, $hi);

    my $aln = $self->build_base_alignment;

    return undef  unless $aln->all_ids > 0;

    if ($self->{'mode'} eq 'new') {
        $aln = $self->build_new_alignment($aln);
    }
    $aln;
}

sub build_indices {
    my $self = shift;
    my ($i, $r, @id);

    $self->{'uid2row'}  = {};
    $self->{'keep_uid'} = {};
    $self->{'hide_uid'} = {};
    $self->{'nops_uid'} = {};

    #index the row objects by unique 'uid' for fast lookup.
    foreach $i (@{$self->{'index2row'}}) {
	$self->{'uid2row'}->{$i->uid} = $i;
    }
    
    #get the reference row handle, if any
    if (@id = $self->map_id($self->{'ref_id'})) {
	$self->{'ref_row'} = $id[0];
    }

    #make all skiplist rows invisible; this has to be done because some
    #may not really have been discarded at all, eg., reference row.
    foreach $i (@{$self->{'skiplist'}}) {
	@id = $self->map_id($i);
	foreach $r (@id) {
	    $self->{'hide_uid'}->{$r->uid} = 1;           #invisible
	}
    }

    #hash the keeplist and make all keeplist rows visible again
    foreach $i (@{$self->{'keeplist'}}) {
	@id = $self->map_id($i);
	foreach $r (@id) {
	    $self->{'keep_uid'}->{$r->uid} = 1;
	    delete $self->{'hide_uid'}->{$r->uid}  if
		exists $self->{'hide_uid'}->{$r->uid};    #visible
	}
    }

    #hash the reference row on the keeplist. don't override
    #any previous invisibility set by discard list.
    $self->{'keep_uid'}->{$self->{'ref_row'}->uid} = 1
	if defined $self->{'ref_row'};
    
    #hash the nopslist: the 'uid' key is used so that the
    #underlying Align class can recognise rows. don't override any previous
    #visibility set by discard list.
    foreach $i (@{$self->{'nopslist'}}) {
	@id = $self->map_id($i);
	foreach $r (@id) {
	    $self->{'nops_uid'}->{$r->uid}  = 1;
	}
    }
    #warn "ref:  ",$self->{'ref_row'}->uid, "\n" if defined $self->{'ref_row'};
    #warn "keep: [", join(",", sort keys %{$self->{'keep_uid'}}), "]\n";
    #warn "nops: [", join(",", sort keys %{$self->{'nops_uid'}}), "]\n";
    #warn "hide: [", join(",", sort keys %{$self->{'hide_uid'}}), "]\n";
    $self;
}

sub build_rows {
    my ($self, $lo, $hi) = @_;

    if ($self->{'aligned'}) {  #treat as alignment: common range
        #warn "range ($lo, $hi)\n";
        for (my $i=0; $i < @{$self->{'index2row'}}; $i++) {
            $self->{'index2row'}->[$i]->assemble($lo, $hi, $self->{'gap'});
        }

    } else {  #treat as format conversion: each row has own range
        for (my $i=0; $i < @{$self->{'index2row'}}; $i++) {
            my ($lo, $hi) = $self->get_range($self->{'index2row'}->[$i]);
            $self->{'index2row'}->[$i]->assemble($lo, $hi, $self->{'gap'});
        }
    }
    $self;
}

sub get_range {
    my ($self, $row) = @_;
    my ($lo, $hi) = $row->range;
    if (@{$self->{'range'}} and @{$self->{'range'}} % 2 < 1) {
	if ($self->{'range'}->[0] < $self->{'range'}->[1]) {
	    ($lo, $hi) = ($self->{'range'}->[0], $self->{'range'}->[1]);
	} else {
	    ($lo, $hi) = ($self->{'range'}->[1], $self->{'range'}->[0]);
	}
    }
    ($lo, $hi);
}

sub build_base_alignment {
    my $self = shift;
    my ($i, $row, $aln, @list) = ();
	
    for ($i=0; $i < @{$self->{'index2row'}}; $i++) {
	$row = $self->{'index2row'}->[$i];
	if (defined $row->{'type'}) {
	    $row = new Bio::MView::Align::Sequence($row->uid, $row->sob,
						   $row->{'type'});
	} else {
	    $row = new Bio::MView::Align::Sequence($row->uid, $row->sob);
	}
	push @list, $row;
    }

    $aln = new Bio::MView::Align(\@list, $self->{'aligned'});

    $aln->set_parameters('nopshash' => $self->{'nops_uid'},
			 'hidehash' => $self->{'hide_uid'});

    #filter alignment based on pairwise %identity, if requested
    if ($self->{'minident'} > 0 or $self->{'maxident'} < 100) {
	$aln = $aln->prune_all_identities($self->{'pcid'},
					  $self->{'minident'},
					  $self->{'maxident'},
					  $self->{'show'},
					  keys %{$self->{'keep_uid'}});
    }

    # foreach my $r ($aln->all_ids) {
    #     $aln->id2row($r)->seqobj->print;
    # }
    # warn "LEN: ", $aln->length;
    $aln;
}

sub build_new_alignment {
    my ($self, $aln) = @_;
    my ($i, $mrow, $arow);

    $aln->set_identity($self->{'ref_row'}->uid, $self->{'pcid'})
	if defined $self->{'ref_row'};

    for ($i=0; $i < @{$self->{'index2row'}}; $i++) {

	$mrow = $self->{'index2row'}->[$i];

	next  if exists $self->{'hide_uid'}->{$mrow->uid};

	$arow = $aln->item($mrow->uid);

	next  unless defined $arow;

	if (exists $self->{'nops_uid'}->{$mrow->uid} or
	    (defined $mrow->{'type'} and $mrow->{'type'} eq 'special')) {
	    $arow->set_display('label0' => '',
			       'label1' => $mrow->cid,
			       'label2' => $mrow->text,
			       'label3' => '',
			       'label4' => '',
			       'label5' => $mrow->posn1,
			       'label6' => $mrow->posn2,
			       'url'    => $mrow->url,
		);
	} else {
	    #don't change label4 (percent identity) here as this was
	    #already computed by $aln->set_identity() above.
	    $arow->set_display('label0' => $mrow->num,
			       'label1' => $mrow->cid,
			       'label2' => $mrow->text,
			       'label3' => $mrow->data,
			       #'label4' => '',###
			       'label5' => $mrow->posn1,
			       'label6' => $mrow->posn2,
			       'url'    => $mrow->url,
		);
	}
    }
    $aln;
}

#remove query and hit columns at gaps in the query sequence and downcase
#the bounding hit symbols in the hit sequence thus affected.
sub strip_query_gaps {
    my ($self, $query, $sbjct) = @_;
    my $i;

    #warn "sqg(in  q)=[$$query]\n";
    #warn "sqg(in  h)=[$$sbjct]\n";

    #no gaps in query
    return    if index($$query, '-') < 0;
    
    #iterate over query frag symbols
    while ( ($i = index($$query, '-')) >= 0 ) {
	
	#downcase preceding symbol in hit
	if (defined substr($$query, $i-1, 1)) {
	    substr($$sbjct, $i-1, 1) = lc substr($$sbjct, $i-1, 1);
	}
	
	#consume gap symbols in query and hit
	while (substr($$query, $i, 1) eq '-') {
	    substr($$query, $i, 1) = "";
	    substr($$sbjct, $i, 1) = "";
	}
	
	#downcase succeeding symbol in hit
	if (defined substr($$query, $i, 1)) {
	    substr($$sbjct, $i, 1) = lc substr($$sbjct, $i, 1);
	}
	
	#warn "sqg(out q)=[$$query]\n";
	#warn "sqg(out h)=[$$sbjct]\n";
    }
    $self;
}


###########################################################################
1;
