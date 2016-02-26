# -*- perl -*-
# Copyright (C) 1996-2015 Nigel P. Brown
# $Id: Universal.pm,v 1.28 2015/06/18 20:26:00 npb Exp $

######################################################################
package Universal;

#useful general stuff
$::Date = `date`;
$::Prog = basename($0);

use strict;

sub member {
    my ($pattern, $list) = @_;
    foreach my $i (@$list) {
        return 1  if $i eq $pattern;
    }
    return 0;
}

#dump sorted list of all instance variables, or supplied list
sub examine {
    my $self = shift;
    my @keys = @_ ? @_ : sort keys %$self;
    print STDERR "Class $self\n";
    foreach my $k (@keys) {
        printf STDERR "%16s => %s\n", $k,
            defined $self->{$k} ? $self->{$k} : '';
    }
    $self;
}

#shallow copy
sub copy {
    my $self = shift;
    my $copy = {};
    foreach my $k (keys %$self) {
	#warn "$k => @{[defined $self->{$k} ? $self->{$k} : 'undef']}\n";
	if (defined $self->{$k}) {
	    $copy->{$k} = $self->{$k};
	} else {
	    $copy->{$k} = '';
	}
    }
    bless $copy, ref $self;
}

#deep copy
sub deep_copy {
    my $self = shift;
    my $copy = {};
    foreach my $k (keys %$self) {
	#warn "$k => @{[defined $self->{$k} ? $self->{$k} : 'undef']}\n";
	if (defined $self->{$k}) {
	    if (my $type = ref $self->{$k}) {
		if (member($type, [ qw(SCALAR ARRAY HASH CODE) ])) {
		    $copy->{$k} = $self->{$k};
		} else {
		    $copy->{$k} = $copy->{$k}->deep_copy;
		}
	    }
	    $copy->{$k} = $self->{$k};
	} else {
	    $copy->{$k} = '';
	}
    }
    bless $copy, ref $self;
}

#warn with error string
sub warn {
    my $self = shift;
    chomp $_[$#_];
    if (ref($self)) {
	warn "Warning ", ref($self), '::', @_, "\n";
	return;
    }
    warn "Warning ", $self, '::', @_, "\n";
}

#exit with error string
sub die {
    my $self = shift;
    chomp $_[$#_];
    if (ref($self)) {
	die "Died ", ref($self), '::', @_, "\n";
    }
    die "Died ", $self, '::', @_, "\n";
}

#replacement for /bin/basename
sub basename {
    my ($path, $ext) = (@_, "");
    if ($^O ne 'MSWin32') {
        ($path) = "/$path" =~ /.*\/(.+)$/;
        return $1  if $path =~ /(.*)$ext$/;
        return $path;
    }
    ($path) = "\\$path" =~ /.*\\(.+)$/;
    return $1  if $path =~ /(.*)$ext$/;
}

#basename and extension
sub fileparts {
    my ($path, $wantbase) = (@_, 1); #discard leading path if true (default)
    $path = basename($path)  if $wantbase;
    return ($1, $2)  if $path =~ /^(.+?)\.([^.]+)$/; #non-greedy
    return ('', $1)  if $path =~ /^\.([^.]+)$/;
    return ($1, '')  if $path =~ /^(.+)\.$/;
    return ($path, '');
}

#temporary file name
sub tmpfile {
    my ($s) = (@_, $$);
    return "/tmp/$s"  if $^O ne 'MSWin32';
    return $s;
}

#arithmetic min() function
sub min {
    my ($a, $b) = @_;
    $a < $b ? $a : $b;
}

#arithmetic max() function
sub max {
    my ($a, $b) = @_;
    $a > $b ? $a : $b;
}

#Linux only?
sub vmstat {
    my ($s) = (@_, '');
    local ($_, *TMP);
    if (open(TMP, "cat /proc/$$/stat|")) {
	$_=<TMP>; my @ps = split /\s+/; close TMP;
	CORE::warn sprintf "VMEM=%8gk  $s\n", $ps[22] / 1024;
    } else {
	CORE::warn sprintf "VMEM=?  $s\n";
    }
}


######################################################################
1;
