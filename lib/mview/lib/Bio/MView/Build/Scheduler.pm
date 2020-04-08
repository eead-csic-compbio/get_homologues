# Copyright (C) 2015-2017 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Build::Scheduler;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::MultiScheduler);

###########################################################################
package Bio::MView::Build::SimpleScheduler;

# An iterator over a stored list of prescribed values (supplied to the
# constructor) which, given an optional subset of these (supplied with the
# filter() method), emits the latter items one by one each time the next()
# method is called, until done when it emits undef.
#
# If the constructor is called with no list of prescribed items, the next()
# method will return a true value just once, effecting a one-pass scheduler.
#
# use Bio::MView::Build::SimpleScheduler;
#
# $s = new Bio::MView::Build::SimpleScheduler;
# $s->next;  #returns each of 1, undef
#
# $s = new Bio::MView::Build::SimpleScheduler([1,2,3]);
# $s->next;  #returns each of 1, 2, 3, undef
# $s->reset;
# $s->next;  #repeats
#
# $s = new Bio::MView::Build::SimpleScheduler([1,2,3]);
# $s->filter([2, 3, 'fred']);
# $s->next;  #returns each of 2, 3, undef

my $DEBUG = 0;

sub new {
    my ($type, $items) = (shift, shift);
    $items = [1]  unless defined $items;  #simple scheduler: returns once
    warn "Scheduler::new:   items [@$items]\n"  if $DEBUG;

    my $self = {}; bless $self, $type;

    $self->{'items'}   = $self->_get_items($items);  #ordered set
    $self->{'todo'}    = undef;                      #values actually wanted
    $self->{'todoidx'} = undef;                      #current todoset index

    $self->filter;

    $self;
}

######################################################################
# public methods
######################################################################
#set wanted values or all of them if no parameter
sub filter {
    my ($self, $todo) = (@_);
    unless (defined $todo) {
        warn "Scheduler::filter: todo undefined - doing all\n"  if $DEBUG > 1;
        $todo = [];  #make an empty list
    }
    $self->{'todo'} = $self->_get_todo($todo);
    warn "Scheduler::filter: todo [@{$self->{todo}}]\n"  if $DEBUG;
    $self->reset;
}

#return 1-based counter
sub itemnum {
    return 0  unless defined $_[0]->{'todoidx'};
    return $_[0]->{'todoidx'} + 1;
}

#return the currently selected item
sub item {
    return undef  unless defined $_[0]->{'todoidx'};
    return $_[0]->{'todo'}->[ $_[0]->{'todoidx'} ];
}

#test if the currently selected item matches the argument
sub use_item {
    return 1  if $_[0]->item eq $_[1];
    return 0;
}

#return the next item and advance the iterator, or return undef after
#completion of a cycle
sub next {
    my $self = shift;
    warn "Scheduler::next entry @{[defined $self->{'todoidx'} ? $self->{'todoidx'} : 'undef']}\n"  if $DEBUG > 1;

    #finished previously: caller needs to reset()
    if ($self->is_done) {
        warn "Scheduler::next finished previously\n"  if $DEBUG;
        return undef;
    }

    #still going - return next value
    if ($self->is_running) {
        $self->{'todoidx'}++;
        my $item = $self->{'todo'}->[ $self->{'todoidx'} ];
        warn "Scheduler::next => $item\n"  if $DEBUG;
        return $item;
    }

    #finished now
    warn "Scheduler::next done\n"  if $DEBUG;
    return $self->{'todoidx'} = undef;
}

#reset iterator back to start
sub reset {
    my $self = shift;
    warn "Scheduler::reset\n"  if $DEBUG;
    $self->{'todoidx'} = -1;  #ready to go
}

#scheduler state: finished waiting for reset
sub is_done {
    return 0  if defined $_[0]->{'todoidx'};
    return 1;
}

#scheduler state: post-reset, ready to start
sub is_ready {
    return 0  unless defined $_[0]->{'todoidx'};
    return 1  if $_[0]->{'todoidx'} == -1;
    return 0;
}

#scheduler state: iterating
sub is_running {
    my $i = $_[0]->{'todoidx'};
    return 0  unless defined $i;
    return 1  if -2 < $i and $i < @{$_[0]->{'todo'}};
}

######################################################################
# private methods
######################################################################
#return list of unique items in given order
sub _get_items {
    my ($self, $items) = @_;
    my ($list, $hash) = ([], {});
    foreach my $i (@$items) {
        next  if exists $hash->{$i};
        push @$list, $i;
        $hash->{$i}++;
    }
    return $list;
}

#return a list of unique and ordered values from $todo as found in
#$self->{'items'}. if $self->{'items'} is in fact empty, return all todo
#elements. the special value 0 in $todo means the 'last' item in
#$self->{'items'}.
sub _get_todo {
    my ($self, $todo) = @_;
    my $range = $self->{'items'};
    warn "Scheduler::_get_todo: entry: items [@$range] todo [@$todo]\n"
        if $DEBUG > 2;

    #empty range: return everything on todo
    if (! @$range and @$todo > 0) {
        my @items = (); my %tmp = ();
        map { $tmp{$_}++ } @$todo;
        foreach my $i (@$todo) {
            push @items, $i  if exists $tmp{$i};
            delete $tmp{$i};
        }
        return \@items;
    }

    #empty todo: return everything on range
    unless (@$todo) {
        warn "Scheduler::_get_todo:  exit: items [@$range] todo [@$range]\n"
            if $DEBUG > 2;
        return [ @$range ];
    }

    #make a dict of requested items for uniqueness
    my @items = (); my %tmp = ();
    map { $tmp{$_}++ } @$todo;
    $tmp{$range->[-1]}++  if exists $tmp{0};  #0 means the 'last' one

    #convert numerical requested items to an indexed value
    foreach my $i (keys %tmp) {
        if ($i =~ /^\d+$/ and 0 < $i and $i <= @$range) {
            delete $tmp{$i};
            $tmp{$range->[$i-1]} = 1;
        }
    }

    #extact legitimate items in correct order
    foreach my $i (@$range) {
        push @items, $i  if exists $tmp{$i};
    }
    warn "Scheduler::_get_todo:  exit: items [@$range] todo [@items]\n"
        if $DEBUG > 2;

    return \@items;
}

###########################################################################
package Bio::MView::Build::MultiScheduler;

# An iterator over a multiple stored lists of prescribed values (supplied to
# the constructor) which, given an optional subset of these (supplied with
# the filter() method), emits the latter items one by one each time the
# next() method is called, until done when it emits undef.

# If the constructor is called with no list of prescribed items, the next()
# method will return a true value just once, effecting a one-pass scheduler.

# use Bio::MView::Build::MultiScheduler;
#
# $s = new Bio::MView::Build::MultiScheduler;
# $s->next;  #returns each of 1, undef
#
# $s = new Bio::MView::Build::MultiScheduler([1,2,3]);
# $s->filter([2, 3, 'fred']);
# $s->next;  #returns each of 2, 3, undef
#
# $s = new Bio::MView::Build::MultiScheduler([1,2,3] ['X','Y']);
# $s->next;  #returns each of [1,X], [1,Y], [2,X], [2,Y] ... [3,Y], undef
# $s->reset;
# $s->next;  #repeats

sub new {
    my $type = shift;
    push @_, [1]  unless @_;
    warn "MultiScheduler::new [@_]\n"  if $DEBUG;

    my $self = {}; bless $self, $type;

    $self->{stage} = [];  #list of scheduler stages

    foreach my $stage (@_) {
        push @{$self->{stage}}, new Bio::MView::Build::SimpleScheduler($stage);
    }
    $self;
}

######################################################################
# public methods
######################################################################
sub filter {
    my $self = shift;
    push @_, []  unless @_;
    die "MultiScheduler::filter: unexpected number of filters\n"
        unless scalar @{$self->{stage}} == scalar @_;
    for (my $i=0; $i<@_; $i++) {
        my $stage = $self->{stage}->[$i];
        $stage->filter($_[$i]);
    }
}

sub itemnum {
    my $self = shift;
    my @tmp = ();
    foreach my $stage (@{$self->{stage}}) {
        push @tmp, $stage->itemnum;
    }
    return @tmp  if wantarray;
    return $tmp[0];  #behave like SimpleScheduler
}

sub item {
    my $self = shift;
    my @tmp = ();
    foreach my $stage (@{$self->{stage}}) {
        push @tmp, $stage->item;
    }
    return @tmp  if wantarray;
    return $tmp[0];  #behave like SimpleScheduler
}

sub use_item {
    my $self = shift;
    die "MultiScheduler::use_item: unexpected number of tests\n"
        unless scalar @{$self->{stage}} == scalar @_;
    for (my $i=0; $i<@_; $i++) {
        my $stage = $self->{stage}->[$i];
        return 0  unless $stage->use_item($_[$i]);
    }
    return 1;
}

#iterate over multiple lists
sub next {
    my $self = shift;
    warn "MultiScheduler::next entry @{[scalar @{$self->{stage}}]} stages\n"
        if $DEBUG;

    for (my $i=@{$self->{stage}}-1; $i>-1; $i--) {

        my $stage = $self->{stage}->[$i];

        warn "  $i: [@{$stage->{todo}}]\n"  if $DEBUG;

        my $r = $stage->next;

        if ($DEBUG) {
            my $num = $stage->itemnum; $num = 'undef' unless defined $num;
            my $itm = $stage->item;    $itm = 'undef' unless defined $itm;
            warn "  $i: $num, $itm\n";
        }

        #this stage yields something
        if (defined $r) {
            #restart any higher stages
            for (my $j=$i-1; $j>-1; $j--) {
                my $stage = $self->{stage}->[$j];
                if ($stage->is_done) {
                    $stage->reset;
                }
                if ($stage->is_ready) {
                    $stage->next;
                }
            }
            return $r;
        }

        #this stage is finished: restart it and go back one, but never
        #reset the first stage, because that marks the end
        if ($i > 0) {
            $stage->reset;
            $stage->next;
            next;
        }
    }
    #all done
    return undef;
}

sub reset {
    my $self = shift;
    warn "MultiScheduler::reset\n"  if $DEBUG;
    foreach my $stage (@{$self->{stage}}) {
        $stage->reset;
    }
}

###########################################################################
1;
