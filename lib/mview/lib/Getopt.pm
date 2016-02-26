# -*- perl -*-
# Copyright (C) 1999-2015 Nigel P. Brown
# $Id: Getopt.pm,v 1.12 2015/06/14 17:09:03 npb Exp $

###########################################################################
# 
# if 'usage' is undefined: don't output any usage information.
# if 'param' is undefined: don't copy option to parameter.
# if 'param' is defined but null, copy with parameter name = option name.
# if 'param' is defined and non-null, copy with new parameter name given
#   by 'param'..
# if 'type'  is empty: treat option as a switch.
# if 'convert' references a function, use any return value from that
#   function to set the parameter. 
# 
###########################################################################
package Getopt::Class;

use Getopt::Long;
use NPB::Parse::Regexps;
use strict;

my $DEBUG = 0;

my %Template =
    (
     'usage'   => undef,
     'type'    => '',
     'default' => undef,
     'param'   => undef,
     'convert' => undef,
     'action'  => undef,
     'print'   => undef,
    );

sub new {
    my $type = shift;
    my ($prog, $options, $passthru) = (@_, 1);
    my $self = {};

    $self->{'prog'}   = $prog;
    $self->{'option'} = $options;
    $self->{'order'}  = [];
    $self->{'errors'} = [];

    #want to search arglist for known options
    Getopt::Long::config(qw(permute));

    if ($passthru) {
	#keep quiet about unknown options: may be
	#recognised by another instance
	Getopt::Long::config(qw(pass_through));
    } else {
	#complain about unknown options
	Getopt::Long::config(qw(no_pass_through));
    }

    bless $self, $type;
}

sub init {
    my $self = shift;
    my ($key, $field);
    foreach $key (keys %{$self->{'option'}->[1]}) {
	foreach $field (keys %Template) {
	    $self->{'option'}->[1]{$key}->{$field} = $Template{$field}
	        if ! exists $self->{'option'}->[1]{$key}->{$field};
	}
    }
    $self;
}

sub set_text {
    my ($self, $text) = @_;
    $self->{'option'}->[0] = $text;
    $self;
}

sub set_option {
    my ($self, $option) = @_;
    $self->{'option'}->[1]->{$option} = {};
    push @{$self->{'order'}}, $option;
    $self;
}

sub set_option_keyval {
    my ($self, $option, $key, $val) = @_;
    #warn "($option, $key, $val)\n";
    $self->{'option'}->[1]->{$option}->{$key} = $val;
    $self;
}

sub warn {
    my $self = shift;
    push @{$self->{'errors'}}, "$self->{'prog'}: @_\n";
    $self;
}

sub die {
    my $self = shift;
    CORE::warn $self->errors    if @{$self->{'errors'}};
    CORE::die "$self->{'prog'}: @_\n";
}

sub usage {
    my $self = shift;
    my ($s, $o, $default, $text, $type) = ('');
    my @list = ();

    push @_, @{$self->{'order'}}    unless @_;
    
    foreach $o (@_) {
	next    unless defined $self->{'option'}->[1]{$o}->{'usage'};
	push @list, $o; 
    }

    if (@list) {
	if (defined $self->{'option'}->[0] and $self->{'option'}->[0]) {
	    $s = $self->{'option'}->[0] . "\n";
	}
    } else {
	return '';
    }
    
    foreach $o (@list) {
	#warn "($o)\n";
	next    unless defined $o;
	next    unless exists $self->{'option'}->[1]{$o};

	$default = $self->{'option'}->[1]{$o}->{'default'};
	if (defined $default) {
	    if ($self->{'option'}->[1]{$o}->{'type'} eq '') {
		$default = ($default ? 'set' : 'unset');
	    } elsif ($default eq '') {
		$default = "''";
	    }
	    $default = "[$default]";
	} else {
	    $default = "[no default]";
	}
	$text = $self->{'option'}->[1]{$o}->{'usage'};
	$type = '';
	$type = "integer"   if $self->{'option'}->[1]{$o}->{'type'} eq 'i';
	$type = "float"     if $self->{'option'}->[1]{$o}->{'type'} eq 'f';
	$type = "string"    if $self->{'option'}->[1]{$o}->{'type'} eq 's';
	$type = "int[,int]" if $self->{'option'}->[1]{$o}->{'type'} eq '@i';
	$type = "flo[,flo]" if $self->{'option'}->[1]{$o}->{'type'} eq '@f';
	$type = "str[,str]" if $self->{'option'}->[1]{$o}->{'type'} eq '@s';
	$type = "int[,int]" if $self->{'option'}->[1]{$o}->{'type'} eq '@I';
	$type = "flo[,flo]" if $self->{'option'}->[1]{$o}->{'type'} eq '@F';
	$type = "str[,str]" if $self->{'option'}->[1]{$o}->{'type'} eq '@S';
	$type = "on|off"    if $self->{'option'}->[1]{$o}->{'type'} eq 'b';
	$s .= sprintf("  -%-20s %s %s.\n", "$o $type", $text, $default);
    }
    return "$s\n";
}

sub option_value {
    my $self = shift;
    my ($o, $v) = @_;
    #warn "($o, $v)\n";
    return "'undef'"    unless defined $v;
    return "''"         if $v eq '';
    return $v;
}

sub param_value {
    my $self = shift;
    my ($o, $v) = @_;
    #warn "($o, $v)\n";
    return "'undef'"    unless defined $v;
    if (ref $v or $self->{'option'}->[1]{$o}->{'type'} =~ /@/) {
	#(assume) it's a list!
	$v = "[" . join(",", @$v) . "]" 
    }
    $v;
}

sub getoptions {
    my $self = shift;
    my ($caller, $opt, $par) = @_;
    my (@tmp, $o, $ov, $p, $pv);

    return    if (@tmp = _build_options($self->{'option'}->[1])) < 1;

    GetOptions($opt, @tmp);

    #map { print STDERR "$_ => $opt->{$_}\n" } %$opt;

OPTION:
    foreach $o (@{$self->{'order'}}) {

	#get raw option value
	if (defined $opt->{$o}) {
	    #from command line
	    $ov = $opt->{$o};
	} else {
	    #use default
	    $ov = $self->{'option'}->[1]{$o}->{'default'};
	}

	#get the parameter name, if any
	$p = $self->{'option'}->[1]{$o}->{'param'};

	#use option name if parameter name was set but empty
	$p = $o    if defined $p and $p eq '';
	
	#perform trivial type tests and trivial option value
	#to parameter conversion
	$pv = $self->test_type($self->{'option'}->[1]{$o}->{'type'}, $o, $ov);

#	if (defined $p) {
#	    warn "option:$o\nparam:$p\novalue:@{[(defined $ov?$ov:'undef')]}\npvalue:@{[(defined $pv?$pv:'undef')]}\n\n";
#	} else {
#		warn "option:$o\novalue:@{[(defined $ov?$ov:'undef')]}\npvalue:@{[(defined $pv?$pv:'undef')]}\n\n";
#	}

	next OPTION    if @{$self->{'errors'}};

	#special option value to parameter conversion, if any
	if (defined $p) {
	    if (defined $self->{'option'}->[1]{$o}->{'convert'} and
		ref $self->{'option'}->[1]{$o}->{'convert'} eq 'CODE') {
		#use conversion function
		$pv = &{$self->{'option'}->[1]{$o}->{'convert'}}
		          ($caller, $self, $o, $ov, $p, $pv);
		#warn "CONV($o): $pv\n"    if defined $pv;
	    }
	}

	next OPTION   if @{$self->{'errors'}};

	#action: perform associated action, if any
	if (defined $self->{'option'}->[1]{$o}->{'action'} and
	    ref $self->{'option'}->[1]{$o}->{'action'} eq 'CODE') {
	    $ov = &{$self->{'option'}->[1]{$o}->{'action'}}
	              ($caller, $self, $o, $ov, $p, $pv);
	    #warn "ACTN($o): $ov\n"    if defined $ov;
	}
	
	next OPTION   if @{$self->{'errors'}};

	#store any converted parameter value
	$opt->{$o} = $ov    if ! exists $opt->{$o};

	#store any converted parameter value
	$par->{$p} = $pv    if defined $p;

	#for debugging
	if ($DEBUG) {
	    $ov = $self->option_value($o, $ov);

	    if (defined $p) {
		$p = sprintf("   %15s => %s", $p, $self->param_value($o, $pv));
	    } else {
		$p = '';
	    }
	    
	    if (defined $opt->{$o}) {
		printf STDERR "found:   %15s => %-10s %s\n", $o, $ov, $p;
	    } else {
		printf STDERR "default: %15s => %-10s %s\n", $o, $ov, $p;
	    }
	}
    }

    @{$self->{'errors'}};
}

sub errors { return @{$_[0]->{'errors'}} }

sub test_type {
    my ($self, $type, $o, $v) = @_;
    return $v    unless defined $type and $type ne '';
    return $v    if $type eq 's';
    if ($type eq 'i') {
	unless ($v =~ /^$RX_Sint$/) {
	    $self->warn("bad argument '$o=$v', want integer");
	} 
	return $v;
    }
    if ($type eq 'f') {
	unless ($v =~ /^$RX_Sreal$/) {
	    $self->warn("bad argument '$o=$v', want float");
	} 
	return $v;
    }
    if ($type eq '@i') {
	return $self->expand_integer_list($o, $v, 0);
    }
    if ($type eq '@I') {
	return $self->expand_integer_list($o, $v, 1);    #expect sorted
    }
    if ($type eq '@f') {
	return $self->expand_float_list($o, $v, 0);
    }
    if ($type eq '@F') {
	return $self->expand_float_list($o, $v, 1);      #expect sorted
    }
    if ($type eq '@s') {
	return $self->expand_list($o, $v, 0);
    }
    if ($type eq '@S') {
	return $self->expand_list($o, $v, 1);            #expect sorted
    }
    if ($type eq 'b') {
	return $self->expand_toggle($o, $v);
    }
    CORE::die "Getopt::Class::test_type() unknown type '$type'\n";
}

sub expand_integer_list {
    my ($self, $o, $v, $sortP) = (@_, 0);
    return []    unless defined $v;
    my @tmp = ();
    local $_;
    #warn "expand_integer_list($o, [$v])\n";
    foreach (split /[,\s]+/, $v) {
	next  unless length($_);
	#warn ">>>[$_]";
        #range M\.\.N or M:N
        if (/^($RX_Sint)(?:\.\.|:)($RX_Sint)$/) {
            if ($2 < $1) {
		if ($sortP) {
		    $self->warn("bad integer list range value '$o=$_'");
		    next;
		} else {
		    push @tmp, $2..$1;
		}
            } else {
		push @tmp, $1..$2;
	    }
            next;
        }
        #non-range
        if (/^($RX_Sint)$/ and ! /\.\./ and ! /:/) {
            push @tmp, $1;
            next;
        }
        $self->warn("bad integer list value '$o=$_'");
        return [];
    }
    #warn "expand_integer_list(@tmp)\n";
    return [ sort @tmp ]    if $sortP;
    return [ @tmp ];
}

sub expand_float_list {
    my ($self, $o, $v, $sortP) = (@_, 0);
    return []    unless defined $v;
    my @tmp = ();
    local $_;
    #warn "expand_float_list($o, [$v])\n";
    foreach (split /[,\s]+/, $v) {
	next  unless length($_);
	#warn ">>>[$_]";
        #non-range
        if (/^($RX_Sreal)$/ and ! /\.\./ and ! /:/) {
            push @tmp, $1;
            next;
        }
        $self->warn("bad float list value '$o=$_'");
        return [];
    }
    #warn "expand_float_list(@tmp)\n";
    return [ sort @tmp ]    if $sortP;
    return [ @tmp ];
}

sub expand_list {
    my ($self, $o, $v, $sortP) = (@_, 0);
    return []    unless defined $v;
    my @tmp = ();
    local $_;
    #warn "expand_list($o, [$v])\n";
    foreach (split /[,\s]+/, $v) {
	next  unless length($_);
	#warn ">>>[$_]";
        #integer range M\.\.N or M:N
        if (/^($RX_Sint)(?:\.\.|:)($RX_Sint)$/) {
            if ($2 < $1) {
		if ($sortP) {
		    $self->warn("bad integer list range value '$o=$_'");
		    next;
		} else {
		    push @tmp, $2..$1;
		}
            } else {
		push @tmp, $1..$2;
	    }
            next;
        }
	#non-range: take whole string
	push @tmp, $_;
    }
    #warn "expand_list(@tmp)\n";
    return [ sort @tmp ]    if $sortP;
    return \@tmp;
}

sub expand_toggle {
    my ($self, $o, $v) = @_;
    return 'off'    unless defined $v;
    if ($v ne 'on' and $v ne 'off' and $v ne '0' and $v ne '1') {
	$self->warn("bad value for '$o=$v' want {on,off} or {0,1}");
    }
    return 1    if $v eq 'on' or $v eq '1';
    return 0;
}

sub _build_options {
    my $opt = shift;
    my @opt = ();
    local $_;
    foreach (keys %$opt) {
	if (!defined $opt->{$_}->{'type'} or $opt->{$_}->{'type'} eq '') {
	    push @opt, $_;
	} else {
	    push @opt, "$_=s";
	}
    }
    #warn "OPT: @opt\n";
    @opt;
}

sub _load_options {
    my ($scope, $prog, $stm) = @_;
    my ($tmp, $class, $name, $option);
    my $text = '';
    my @class = ();
    my %class;
    local $_;

    while (<$stm>) {
	
	last    if eof;

	#blank
	next    if /^\s*$/;

	#hash comment
	next    if /^\s*\#/;

	chomp;

	#TEXT
	if (!defined $class and /^\s*text\s*:\s*(.*)/i) {
	    #warn "#text($1)\n";
	    ($text, $_) = _scan_quoted_text($prog, $stm, $1);
	    redo;
	}
	
	#CLASS
	if (/^\s*\[\s*([a-z0-9]+)\s*\]/i) {
	    $name = $1;
	    if (! exists $class{$name}) {
		$class = new Getopt::Class($prog, [], 1);
		$class{$name} = $class;
		push @class, $name;
	    } else {
		CORE::die "Getopt::Class::_load_options() class duplication for '$name'\n";
	    }
	    next;
	}

	#class.TEXT
	if (/^\s*text\s*:\s*(.*)/i) {
	    #warn "#class.text($1)\n";
	    ($tmp, $_) = _scan_quoted_text($prog, $stm, $1);
	    $class->set_text($tmp);
	    redo;
	}
	
	#class.OPTION
	if (/^\s*option\s*:\s*(\S+)/i) {
	    #warn "#class.option($1)\n";
	    $option = _strip_quotes($1);
            $class->set_option($option);
	    next;
	}
	
	#class.option.TYPE
	if (/^\s*(type)\s*:\s*(\S+)/i) {
	    #warn "#class.option.$1($2)\n";
            $class->set_option_keyval($option, $1, _strip_quotes($2));
	    next;
	}
	
	#class.option.DEFAULT
	if (/^\s*(default)\s*:\s*(.*)/i) {
	    #warn "#class.option.$1($2)\n";
	    $class->set_option_keyval($option, $1, _strip_quotes($2));
	    next;
	}
	
	#class.option.USAGE
	if (/^\s*(usage)\s*:\s*(.*)/i) {
	    #warn "#class.option.$1($2)\n";
	    ($tmp, $_) = _scan_quoted_text($prog, $stm, $2);
	    $class->set_option_keyval($option, $1, $tmp);
	    redo;
	}
	
	#class.option.PARAM
	if (/^\s*(param)\s*:\s*(\S*)/i) {
	    #warn "#class.option.$1($2)\n";
	    $class->set_option_keyval($option, $1, _strip_quotes($2));
	    next;
	}
	
	#class.option.CONVERT
	if (/^\s*(convert)\s*:\s*(.*)/i) {
	    #warn "#class.option.$1($2)\n";
	    ($tmp, $_) = _scan_subroutine($scope, $prog, $stm, "$name.$option", $2);
	    $class->set_option_keyval($option, $1, $tmp);
	    redo;
	}
	
	#class.option.ACTION
	if (/^\s*(action)\s*:\s*(.*)/i) {
	    #warn "#class.option.$1($2)\n";
	    ($tmp, $_) = _scan_subroutine($scope, $prog, $stm, "$name.$option", $2);
	    $class->set_option_keyval($option, $1, $tmp);
	    redo;
	}
	
	#class.option.PRINT
	if (/^\s*(print)\s*:\s*(.*)/i) {
	    #warn "#class.option.$1($2)\n";
	    ($tmp, $_) = _scan_subroutine($scope, $prog, $stm, "$name.$option", $2);
	    $class->set_option_keyval($option, $1, $tmp);
	    redo;
	}

	CORE::die "Getopt::Class::_load_options() unrecognised line: [$_]";
    }
    ($text, \@class, \%class);
}

sub _scan_quoted_text {
    my ($prog, $stm, $line, $text) = (shift, shift, shift, '');
    $line = ''    unless defined $line;
    #warn "($stm, $line, $text)";
    if ($line =~ /([\"\'].*[\"\'])\s*$/) {
	$text = $1;                                     #single line
	$line = <$stm>;
    }
    elsif ($line =~ /([\"\'].*)/) {
	$text = $1;                                     #first line
	while ($line = <$stm>) {
	    last    if $line =~ /^\s*\S+\s*:/;          #next option
	    last    if $line =~ /^\s*\[\s*[a-z0-9]+\s*\]i/; #next class
	    chomp $line;
	    if ($line =~ /(.*[\"\'])\s*$/) {            #last line
		$text .= $1;
		$line = <$stm>;
		last;
	    }
	    $text .= $line;                             #middle lines
	}
    }
    #warn "TXT: ($text)\n";
    $text =~ s/^[\"\']//;       #strip leading quote
    $text =~ s/[\"\']$//;       #strip trailing quote
    $text =~ s/\\n/\n/g;        #translate newlines
    $text = _process_macros($prog, $text, 1);
    #warn "TXT: ($text)\n"; 
    ($text, $line);
}

sub _scan_subroutine {
    my ($scope, $prog, $stm, $option, $line) = @_;
    my $tmp = '';
    $line = ''    unless defined $line;
    #warn "($stm, $line, $tmp)";
    if ($line =~ /^\s*(sub.*)/) {
        $tmp = "$1\n";                                  #first line
    }
    while ($line = <$stm>) {    
        last if $line =~ /^\s*(?:option|usage|type|default|param|convert|action|print)\s*:/i;                                          #next option
        last    if $line =~ /^\s*\[\s*[a-z0-9]+\s*\]/i; #next class
        $tmp .= $line;                                  #middle lines
    }

    #warn "SUB: ($tmp)\n";
    $tmp = _process_macros($prog, $tmp, 0);
    #warn "SUB: ($tmp)\n";

    $tmp = eval $tmp;
    CORE::die "Getopt::Class::_load_options() bad sub definition '$option':\n$@"    if $@;

    ($tmp, $line);
}

sub _strip_quotes {
    my $txt = shift;
    $txt =~ s/^[\"\']//;
    $txt =~ s/[\"\']$//;
    $txt;
}

sub _process_macros {
    my ($prog, $text, $string) = @_;

    #__PROG macro
    $text =~ s/__PROG/$prog/g;

    if ($string) {
	#__CHOOSE() macro
	if ($text =~ /__CHOOSE\((.*)\)/) {
	    my ($repl, $sig, @tmp);
	    $sig = $SIG{'__WARN__'};
	    $SIG{'__WARN__'} = sub {};
	    @tmp = eval $1;
	    if ($@) {
		$repl = $1;
	    } else {
		$repl = join(",", @tmp);
	    }
	    $SIG{'__WARN__'} = $sig;
	    $text =~ s/__CHOOSE\((.*)\)/\{$repl\}/g;
	}
    } else {
	#__ARGS macro
	$text =~ s/__ARGS/my (\$getopt,\$class,\$on,\$ov,\$pn,\$pv)=\@_;/g;

	#__GETOPT macro
	$text =~ s/__GETOPT/\$getopt/g;

	#__CLASS macro
	$text =~ s/__CLASS/\$class/g;

	#__OPTION macro
	$text =~ s/__OPTION/\$getopt->{'option'}->/g;
	
	#__PARAM macro
	$text =~ s/__PARAM/\$getopt->{'param'}->/g;
	
	#__DELETE_OPTION macro
	$text =~ s/__DELETE_OPTION/\$getopt->delete_option/g;
	
	#__DELETE_PARAM macro
	$text =~ s/__DELETE_PARAM/\$getopt->delete_parameter/g;
	
	#__TEST macro
	$text =~ s/__TEST/\$class->test_type/g;

	#__WARN macro
	$text =~ s/__WARN/\$class->warn/g;

	#__USAGE macro
	$text =~ s/__USAGE/\$getopt->usage/g;
	
	#__ON macro (option name)
	$text =~ s/__ON/\$on/g;

	#__OV macro (option value)
	$text =~ s/__OV/\$ov/g;

	#__PN macro (parameter name)
	$text =~ s/__PN/\$pn/g;

	#__PV macro (parameter value)
	$text =~ s/__PV/\$pv/g;
    }
    $text;
}


###########################################################################
package Getopt;

use strict;

sub new {
    my $type = shift;
    my ($prog, $stm) = @_;
    my $self = {};
    local $_;

    (
     $self->{'text'},
     $self->{'order'},
     $self->{'class'},
    ) = Getopt::Class::_load_options((caller)[0], $prog, $stm);
    $self->{'prog'}   = $prog;
    $self->{'option'} = {};
    $self->{'param'}  = {};
    $self->{'argv'}   = [];
	
    foreach (keys %{$self->{'class'}}) {
	$self->{'class'}->{$_}->init;
    }

    bless $self, $type;
}

sub usage {
    my ($self, $stm) = (@_, \*STDERR);
    my ($c, $k);
    print $stm "$self->{'text'}\n"    if defined $self->{'text'};
    foreach $c (@{$self->{'order'}}) {
	#print $stm "[$c]\n", $self->{'class'}->{$c}->usage;
	print $stm $self->{'class'}->{$c}->usage;
    }
    $self;
}

sub getoptions {
    my ($self, $argv, $stm) = (@_, \*STDERR);
    my @tmp = ();
    my $error = 0;
    local $_;

    #save input ARGV for posterity
    push @{$self->{'argv'}}, @$argv;

    foreach $_ (@{$self->{'order'}}) {
	if ($self->{'class'}->{$_}->getoptions($self, $self->{'option'}, 
					       $self->{'param'})) {
	    print $stm $self->{'class'}->{$_}->errors;
	    $error++;
	}
    }
    #warn "remaining1 ARGV: @ARGV\n";
    foreach (@ARGV) {
	if (/^--?\S/) {
	    CORE::warn "$self->{'prog'}: unknown or badly formed option '$_'\n";
	    $error++;
	} else {
	    push @tmp, $_;
	}
    }
    CORE::die "$self->{'prog'}: aborting.\n"  if $error;
    @ARGV = @tmp;
    #warn "remaining2 ARGV: @ARGV\n";
    $self;
}

sub get_option_hash { $_[0]->{'option'} }

sub delete_option {
    my $self = shift;
    local $_;
    foreach (@_) {
	delete $self->{'option'}->{$_}    if exists $self->{'option'}->{$_};
    }
    $self;
}

sub argv { @{$_[0]->{'argv'}} }
sub argv_string { join(" ", @{$_[0]->{'argv'}}) }

sub print_options {
    my ($self, $stm) = (@_, \*STDERR);
    local $_;
    print $stm "Option list:\n";
    map { print $stm "$_ => @{[defined $self->{'option'}{$_} ?
        $self->{'option'}{$_} : 'undef']}\n" } sort keys %{$self->{'option'}};
    $self;
}

sub get_parameter_hash { $_[0]->{'param'} }

sub delete_parameter {
    my $self = shift;
    local $_;
    foreach (@_) {
	delete $self->{'param'}->{$_}    if exists $self->{'param'}->{$_};
    }
    $self;
}

sub print_parameters {
    my ($self, $stm) = (@_, \*STDERR);
    local $_;
    print $stm "Parameter list:\n";
    map { print $stm "$_ => @{[defined $self->{'param'}{$_} ?
        $self->{'param'}{$_} : 'undef']}\n" } sort keys %{$self->{'param'}};
    $self;
}


###########################################################################
1;
