# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Parse::Record;

use Bio::Parse::Block;
use Bio::Parse::BlockKeeper;
use Bio::Parse::Scanner;
use Bio::Util::Object;

use vars qw(@ISA);

@ISA = qw(Bio::Util::Object);

my $KEY_DELIM  = "::";

my $IGNORE_ATTR = "parent|text|start|stop|blockkeeper|index|relative_key|absolute_key";

my $DEBUG = 1;

sub new {
    my $ignore = shift;  #discard system supplied type
    my ($type, $parent, $text, $start, $stop) = (@_, -1, -1);

    #warn "> Record($type, $parent, $text, $start, $stop)\n"  if $DEBUG;

    my $self = {};
    bless $self, $type;

    #warn $text->substr(), "\n";

    $self->{'parent'} = $parent;  #parent record
    $self->{'text'}   = $text;    #text object being parsed
    $self->{'start'}  = $start;   #text start
    $self->{'stop'}   = $stop;    #text stop

    $self->{'blockkeeper'} = new Bio::Parse::BlockKeeper();

    #subrecord counter
    $self->{'index'} = $self->get_next_record_number($parent);

    #relative key for reporting
    $self->{'relative_key'} = $self->make_relative_key();

    #absolute hierarchical key for indexing/reporting
    $self->{'absolute_key'} = $self->make_absolute_key();

    #print $self->examine;

    $self;
}

# explicitly call this when finished with a Record to break circular
# references through upward 'parent' and downward 'blockkeeper'.
sub free {
    my $self = shift;

    #warn "> Record::free $self\n"  if $DEBUG;

    #free all subrecords
    $self->{'blockkeeper'}->free();
    $self->{'blockkeeper'} = undef;

    $self->{'text'} = undef;
    $self->{'parent'} = undef;  #break circle

    #warn "< Record::free $self\n"  if $DEBUG;
}

#free records by key
sub free_parsers {
    my $self = shift;
    #warn "> Record::free_parsers $self\n"  if $DEBUG;

    $self->{'blockkeeper'}->free_parsers(@_);

    #warn "< Record::free_parsers $self\n"  if $DEBUG;
}

sub get_offset { $_[0]->{'start'}; }

sub get_bytes { $_[0]->{'stop'} - $_[0]->{'start'}; }

sub get_text { $_[0]->{'text'}; }

# save a subrecord
sub push_record {
    my ($self, $key, $start, $stop) = @_;
    #warn "push_record: incoming=($key, $start, $stop)\n"  if $DEBUG;

    my $block = new Bio::Parse::Block($key, $start, $stop);

    $self->{'blockkeeper'}->add($block);
}

# remove the most recently stored block; returns the deleted block
sub pop_record {
    return $_[0]->{'blockkeeper'}->remove_last();
}

sub get_parent {
    my ($self, $depth) = (@_, 1);
    $depth = -$depth  if $depth < 0;  #positive
    return $self  if $depth == 0;
    my $ptr = $self;
    while ($depth-- > 0) {
        return undef  unless defined $ptr->{'parent'};
        return $ptr->{'parent'}  if $depth < 1;
        $ptr = $ptr->{'parent'};
    }
    return undef;
}

sub get_record {
    my ($self, $key, $num) = (@_, 1);
    my $block = $self->{'blockkeeper'}->get_block($key, $num);
    $self->die("block $key/$num not found")  unless defined $block;
    return $block->parse($self, $self->{'text'});
}

# return list of attributes of Record, excepting housekeeping ones.
sub list_attrs {
    my $self = shift;
    my ($key, @attr);
    foreach $key (grep !/^($IGNORE_ATTR)$/o , keys %$self) {
        push @attr, $key;
    }
    return @attr;
}

sub test_args {
    my ($self, $line) = (shift, shift);
    my $bad = 0; map { $bad++  if $_ eq '' } @_;
    return  unless $bad;
    $line = $$line; chomp $line;
    $self->warn("incomplete match of line: '$line'");
}

sub test_records {
    my $self = shift;
    foreach my $key (@_) {
        $self->warn("corrupt or missing '$key' record\n")
            unless $self->{'blockkeeper'}->has_key($key);
    }
}

# parse all subrecords if called with no argument or '*';
# parse all subrecords of type given by explicit key;
# parse single subrecord of type given by explicit key and 1-based index.
sub parse {
    my ($self, $key, $num) = @_;
    my @list;
    #warn "PARSE(key=@{[defined $key?$key:'']}, num=@{[defined $num?$num:'']}) @{[wantarray]}\n";

    if (@_ == 0 or $key eq '*') {
        my @keys = $self->{'blockkeeper'}->keys();
        @list = $self->parse_many(@keys);
    } elsif (! defined $num) {
        @list = $self->parse_many($key);
    } else {
        @list = $self->parse_one($key, $num);
    }

    #warn "[@list] @{[wantarray]}\n";

    return @list    if wantarray;
    return $list[0] if @list;
    return undef;  #no data
}

sub parse_one {
    my ($self, $key, $num) = (@_, 1);
    #warn " parse_one: $key $num\n";

    my $block = $self->{'blockkeeper'}->get_block($key, $num);
    return ()  unless defined $block;
    return ($block->parse($self, $self->{'text'}));
}

sub parse_many {
    my $self = shift;
    my @list = ();
    foreach my $key (@_) {
        #warn " parse_many: $key\n";
        my @blocks = $self->{'blockkeeper'}->get_block_list($key);
        foreach my $block (@blocks) {
            next  unless defined $block;
            push @list, $block->parse($self, $self->{'text'});
        }
    }
    return @list;
}

# given a record key for this entry, count how many corresponding records
# exist and return an array of counts or the first count if called in a scalar
# context. If called with no argument or just '*' count everything.
sub count {
    my $self = shift;
    my ($key) = @_;
    my (@list, @keys) = ();

    #warn "count(key=@{[defined $key?$key:'']})\n";

    if (@_ == 0 or $key eq '*') {
        @keys = sort $self->{'blockkeeper'}->keys();  #note: sorted
    } else {
        push @keys, $key;
    }

    foreach $key (@keys) {
        push @list, $self->{'blockkeeper'}->count($key);
    }

    return @list    if wantarray;
    return $list[0] if @list;
    return 0;  #no data
}

###########################################################################
# private methods
###########################################################################
#sub DESTROY { warn "DESTROY $_[0]\n" }

sub make_relative_key {
    return $_[0]->basename_class() . $KEY_DELIM . $_[0]->{'index'};
}

sub make_absolute_key {
    my $self = shift;
    my $key = '';
    $key = $self->{'parent'}->{'absolute_key'}. $KEY_DELIM
        if defined $self->{'parent'};
    $key .= $self->{'relative_key'};
    return $key;
}

sub get_next_record_number {
    my ($self, $parent) = @_;
    return 0  unless defined $parent;
    my $i = $parent->{'blockkeeper'}->lookup_block_number($self->basename_class(),
                                                          $self->{'start'});
    return 1  if $i == -1;
    return 1 + $i;
}

# return an array of records indexed by key and optional index (1-based);
# overrides Bio::Util::Object::make_message_string used by
# Bio::Util::Object::warn, Bio::Util::Object::die
sub make_message_string {
    my ($self, $prefix) = (shift, shift);
    my $s = $prefix;
    if (ref $self) {
        my $type = "$self"; $type =~ s/=.*//;
        my $path = $self->{'absolute_key'};
        $s .= " $type ($path)";
    } else {
        $s .= " $self";
    }
    $s .= ": " . Bio::Util::Object::_args_as_string(@_)  if @_;
    return $s;
}

###########################################################################
# debug methods
###########################################################################
# called in the parser test suite
sub print {
    my ($self, $indent) = (@_, 0);
    print $self->dump_record($indent);
}

sub dump_record {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "%sClass:  %s\n", $x, $self;
    $s .= sprintf "%sParent: %s\n", $x,
        defined $self->{'parent'} ? $self->{'parent'} : 'undef';
    $s .= sprintf "%sKey:    %s   Indices: []\n", $x, $self->{'relative_key'};

    $s .= sprintf "%s  Subrecords by posn:\n", $x;
    $s .= $self->{'blockkeeper'}->dump_records_by_posn($indent);

    $s .= sprintf "%s  Miscellaneous:\n", $x;
    $s .= sprintf "$x%20s -> %d\n",   'index', $self->{'index'};
    $s .= sprintf "$x%20s -> [%s]\n", 'pos',
        join(', ', $self->{'start'}, ($self->{'stop'}-$self->{'start'}));
    my $tmp = '';
    if (defined $self->{'text'}) {
        $tmp   = $self->{'text'}->substr($self->{'start'}, 30);
        ($tmp) = split("\n", $tmp);
    }
    $s .= sprintf "$x%20s -> \"%s\" ...\n",  'text', $tmp;

    $s .= sprintf "%s  Data:\n", $x;
    $s .= $self->dump_data($indent);  #supplied by subclass

    return $s;
}

# helper for print(); subclass overrides to add fields
sub dump_data { '' }

# helper for dump_data(); used by subclasses
sub fmt {
    my ($self, $val, $undef) = (@_, '<undef>');
    my $ref = ref $val;
    if ($ref eq 'HASH') {
        my @tmp = map {
            my $v = defined $val->{$_} ? $val->{$_} : $undef; "$_:$v"
        } sort keys %$val;
        return "{" . join(',', @tmp) . "}";
    }
    if ($ref eq 'ARRAY') {
        my @tmp = map { defined $_ ? $_ : $undef } @$val;
        return "[" . join(',', @tmp) . "]";
    }
    return defined $val ? $val : $undef;
}

# given a record key and optional 1-based index, return the subrecord string;
# if called with no argument or '*' return the whole record string.
sub string {
    my $self = shift;
    my ($key, $num) = (@_, 1);

    #warn "string(key=@{[defined $key?$key:'']})\n";

    #whole record string
    if (@_ == 0 or $key eq '*') {
        return $self->{'text'}->substr($self->{'start'},
                                       $self->{'stop'} - $self->{'start'});
    }

    my $block = $self->{'blockkeeper'}->get_block($key, $num);

    return ''  unless defined $block;

    return $self->{'text'}->substr($block->{'start'},
                                   $block->{'stop'} - $block->{'start'});
}

###########################################################################
1;
