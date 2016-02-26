# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Row.pm,v 1.17 2015/06/20 21:55:23 npb Exp $

###########################################################################
package Bio::MView::Build::Row;

use Bio::MView::Sequence;

use strict;

my $DEF_TEXTWIDTH = 30;  #default width to truncate 'text' field
my $DEF_PAD = '-';       #default terminal gap character
my $DEF_GAP = '-';       #default internal gap character

sub new {
    my $type = shift;
    my ($num, $id, $desc) = splice @_, 0, 3;
    my $self = {};

    bless $self, $type;

    #strip non-identifier leading rubbish:  >  or /:
    $id =~ s/^(>|\/:)//;

    $self->{'rid'}  = $id;                      #supplied identifier
    $self->{'uid'}  = $self->uniqid($num, $id); #unique compound identifier
    
    #ensure identifier is non-null (for Build::map_id())
    $id = ' '  unless $id =~ /./;

    #set row 'subtype' information
    if ($id =~ /^\#/) {
	$self->{'type'} = 'special';  #leading hash: user supplied row
	$id =~ s/^\#//;               #and strip it
    } else {
	$self->{'type'} = undef;                    
    }

    $self->{'num'}  = $num;                     #row number/string
    $self->{'cid'}  = $id;                      #cleaned identifier
    $self->{'desc'} = $desc;                    #description string
    $self->{'frag'} = [];                       #list of fragments

    $self->{'seq'}  = new Bio::MView::Sequence; #finished sequence

    $self->{'url'}  = Bio::SRS::srsLink($self->{'cid'});  #url

    $self->{'data'} = {};                       #other parsed info

    $self->save_info(@_);

    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

sub uniqid { "$_[1]\034/$_[2]" }

#methods returning standard strings for use in generic output modes
sub rid  { $_[0]->{'rid'} }
sub uid  { $_[0]->{'uid'} }
sub cid  { $_[0]->{'cid'} }
sub num  { $_[0]->{'num'} }
sub url  { $_[0]->{'url'} }
sub sob  { $_[0]->{'seq'} }

sub desc { $_[0]->{'desc'} }  #row description

sub seq {                    #the sequence
    my ($self, $pad, $gap) = (@_, $DEF_PAD, $DEF_GAP);
    return ''  unless defined $self->{'seq'};
    $self->set_pad($pad);
    $self->set_gap($gap);
    return $self->{'seq'}->string
}

sub posn1 { '' }              #first sequence range
sub posn2 { '' }              #second sequence range

sub text {                    #possibly truncated description
    my $w = defined $_[1] ? $_[1] : $DEF_TEXTWIDTH;
    $w = length $_[0]->{'desc'}  if $w > length $_[0]->{'desc'};
    sprintf("%-${w}s", $_[0]->truncate($_[0]->{'desc'}, $w));
}

#convert nucleotide positions to a relative amino acid scale
sub translate_range {
    my ($self, $fm, $to) = @_;
    return (int(($fm+2)/3), int($to/3))   if $fm < $to;  #orientation +
    return (int($fm/3),  int(($to+2)/3))  if $fm > $to;  #orientation -
    die "translate_range: from == to  $fm, $to";
}

#truncate a string
sub truncate {
    my ($self, $s, $n, $t) = (@_, $DEF_TEXTWIDTH);
    $t = substr($s, 0, $n);
    substr($t, -3, 3) = '...'    if length $s > $n;
    $t;
}

sub print {
    sub _format {
	my ($self, $k, $v) = @_;
	$v = 'undef' unless defined $v;
	$v = "'$v'" if $v =~ /^\s*$/;
	return sprintf("  %-15s => %s\n", $k, $v)
    }
    my $self = shift;
    warn "$self\n";
    map { warn $self->_format($_, $self->{$_}) } sort keys %{$self};
    $self;
}

#routine to sort 'frag' list: default is null
sub sort {$_[0]}

#modify the extra 'type' information
sub set_subtype { $_[0]->{'type'} = $_[1] }

#add a sequence fragment to the 'frag' list with value and positions given
#by first three args. use default positions if called with one arg. other
#optional arguments are special to any subclass of Row.
sub add_frag {
    my $self = shift;
    my ($frag, $qry_from, $qry_to) = (shift, shift, shift);

    $qry_from = 1               unless defined $qry_from;
    $qry_to   = length $frag    unless defined $qry_to;

    push @{$self->{'frag'}}, [ \$frag, $qry_from, $qry_to, @_ ];

    #warn "@{$self->{'frag'}->[-1]}\n";

    $self;
}

sub count_frag { scalar @{$_[0]->{'frag'}} }

#compute the maximal positional range of a row
sub range {
    my $self = shift;
    my ($lo, $hi) = ($self->{'frag'}->[0][1], $self->{'frag'}->[0][2]);
    foreach my $frag (@{$self->{'frag'}}) {
        #warn "range: $frag->[1], $frag->[2]\n";
        $lo = $frag->[1]  if $frag->[1] < $lo;
        $lo = $frag->[2]  if $frag->[2] < $lo;
	$hi = $frag->[1]  if $frag->[1] > $hi;
	$hi = $frag->[2]  if $frag->[2] > $hi;
    }
    #warn "range: ($lo, $hi)\n";
    ($lo, $hi);
}

#assemble a row from sequence fragments
sub assemble {
    my ($self, $lo, $hi, $gap) = @_;
    my $reverse = 0;
    #get direction from first fragment range longer than 1
    foreach my $frag (@{$self->{'frag'}}) {
        $reverse = 0, last  if $frag->[1] < $frag->[2];
        $reverse = 1, last  if $frag->[1] > $frag->[2];
    }
    #warn "Row::assemble: [@_] $reverse\n";
    $self->sort;                                 #fragment order
    $self->{'seq'}->reverse  if $reverse;        #before calling insert()
    $self->{'seq'}->insert(@{$self->{'frag'}});  #assemble fragments
    $self->{'seq'}->set_range($lo, $hi);         #set sequence range
    $self->{'seq'}->set_pad($gap);
    $self->{'seq'}->set_gap($gap);
    $self;
}

sub set_pad { $_[0]->{'seq'}->set_pad($_[1]) }
sub set_gap { $_[0]->{'seq'}->set_gap($_[1]) }
sub set_spc { $_[0]->{'seq'}->set_spc($_[1]) }

###########################################################################

sub schema { die "@{[ref $_[0]]}: no schema found\n" }

#save row information following a schema
sub save_info {
    my $self = shift;
    #warn "save_info: [@_]\n";

    my $schema = eval { $self->schema };
    return $self  unless defined $schema;

    for (my $i=0; $i<@$schema; $i++) {
        my ($n1, $n2, $name, $string, $format, $default) = @{$schema->[$i]};
        #warn "save: $name\n";
        $self->{'data'}->{$name} = defined $_[$i] ? $_[$i] : $default;
    }
    $self;
}

#set a row information attribute if in the schema
sub set_val {
    my ($self, $key, $val) = @_;
    #warn "set_val: [$key => $val]\n";

    my $schema = eval { $self->schema };
    return $self  unless defined $schema;

    for (my $i=0; $i<@$schema; $i++) {
        my ($n1, $n2, $name, $string, $format, $default) = @{$schema->[$i]};
        $self->{'data'}->{$name} = $val, return $self  if $key eq $name;
        $self->{'data'}->{$name} = $val, return $self  if $key eq $string;
    }
    warn "@{[ref $self]}::set_val: unknown attribute '$key'\n";
    $self;
}

#test if item has attribute
sub has {
    my ($self, $key) = @_;
    #warn "has: [$key]\n";

    my $schema = eval { $self->schema };

    if (! defined $schema) {
        return exists $self->{$key};
    }

    for (my $i=0; $i<@$schema; $i++) {
        my ($n1, $n2, $name, $string, $format, $default) = @{$schema->[$i]};
        return 1  if $key eq $name;
        return 1  if $key eq $string;
    }
    0; #no key
}

#get a row information attribute if in the schema
sub get_val {
    my ($self, $key) = @_;
    #warn "get_val: [$key]\n";

    my $schema = eval { $self->schema };

    if (! defined $schema) {
        return $self->{$key}  if exists $self->{$key};
        return '';
    }

    for (my $i=0; $i<@$schema; $i++) {
        my ($n1, $n2, $name, $string, $format, $default) = @{$schema->[$i]};
        return $self->{'data'}->{$name}  if $key eq $name;
        return $self->{'data'}->{$name}  if $key eq $string;
    }
    warn "@{[ref $self]}::get_val: unknown attribute '$key'\n";
    '';
}

#return a string of formatted row information following the schema
sub data {
    my ($self, $delim) = (@_, ' ');

    my $fmtstr = sub {
        my $fmt = shift;
        $fmt =~ /(\d+)(\S)/o;
        "%$1s";
    };

    my $schema = eval { $self->schema };
    return ''  unless defined $schema;

    my @tmp = ();
    for (my $i=0; $i<@$schema; $i++) {
        my ($n1, $n2, $name, $string, $format, $default) = @{$schema->[$i]};
        next  unless $n1 > 0;  #ignore this datum
        my $fmt = &$fmtstr($format);
        #warn "data: $name\n";
        my $data = $self->{'data'}->{$name};
        push(@tmp, sprintf($fmt, $data)),  next  if $self->num;  #data
        push(@tmp, sprintf($fmt, $string))                       #header
    }
    return join($delim, @tmp);
}

#return a list of row information following the schema
sub rdb_info {
    my ($self, $mode) = @_;

    my $schema = eval { $self->schema };
    return ''  unless defined $schema;

    my @tmp = ();
    for (my $i=0; $i<@$schema; $i++) {
        my ($n1, $n2, $name, $string, $format, $default) = @{$schema->[$i]};
        next  unless $n2 > 0;  #ignore this datum
        my $data = $self->{'data'}->{$name};
        #warn "rdb: $name\n";
        push(@tmp, $string), next  if $mode eq 'attr';
        push(@tmp, $format), next  if $mode eq 'form';
        push(@tmp, $data),   next  if $mode eq 'data';
    }
    @tmp;
}

#column labels used in headers and rulers: wrapper for 'sub data'
sub head {''}           #row info labels
sub pcid { '' }         #identity% label
sub pcid_std { 'id%' }  #standard text for identity% label

#tabulated output: called from Convert
sub rdb_row {
    my ($self, $mode, $pad, $gap) = @_;
    my @cols;
    @cols = ('row', 'id', 'desc', 'seq')   if $mode eq 'attr';
    @cols = ('4N', '30S', '500S', '500S')  if $mode eq 'form';
    @cols = ($self->num, $self->cid, $self->desc, $self->seq($pad, $gap))
	if $mode eq 'data';

    my @new = $self->rdb_info($mode);     #subtype has any data?
    splice(@cols, -1, 0, @new)  if @new;  #insert penultimately

    #warn "[@cols]";
    return join("\t", @cols);
}


###########################################################################
package Bio::MView::Build::Simple_Row;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row);

sub new {
    my $type = shift;
    my ($num, $id, $desc, $seq) = @_;
    my $self = new Bio::MView::Build::Row($num, $id, $desc);
    bless $self, $type;
    $self->add_frag($seq)  if defined $seq;
    $self;
}


###########################################################################
1;
