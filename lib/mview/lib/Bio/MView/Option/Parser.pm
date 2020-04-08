# Copyright (C) 2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Option::Parser;

use Getopt::Long qw(GetOptionsFromArray);
use Bio::MView::Option::Parameters;

my @OptionTypeLibs = ( 'Bio::MView::Option::Types' );
my @OptionLibs     = ( 'Bio::MView::Option::Options' );

my $DEBUG=0;

###########################################################################
sub new {
    my ($type, $prog) = @_;
    my $self = {};
    bless $self, $type;

    $self->{'param'}   = { 'prog' => $prog, 'argv' => [] };

    $self->{'types'}   = {};
    $self->{'options'} = {};
    $self->{'groups'}  = [];

    foreach my $lib (@OptionTypeLibs) {
        load_library($lib);
        $self->load_types($lib);
    }

    my $groups = {};
    foreach my $lib (@OptionLibs) {
        load_library($lib, $groups);
        $self->load_option_groups($lib);
    }

    if (! defined $Bio::MView::Option::Options::Header) {
        die "Options: no options header found\n";
    }

    #warn "groups:  ", scalar @{$self->{'groups'}}, "\n";
    #warn "options: ", scalar keys %{$self->{'options'}}, "\n";
    #warn "types:   ", scalar keys %{$self->{'types'}}, "\n";

    #want to search arglist for known options
    Getopt::Long::config(qw(permute));

    #keep quiet about unknown options
    Getopt::Long::config(qw(pass_through));

    return $self;
}

######################################################################
# public methods
######################################################################
sub parse_argv {
    my ($self, $argv) = @_;
    my @errors = ();

    #save input ARGV for posterity
    $self->{'param'}->{'argv'} = [@$argv];

    #look for '--' to stop option processing
    #warn "argv: [@$argv]\n";
    my @lhs = (); my @rhs = (); my $options = 1;
    while (@$argv) {
        my $arg = shift @$argv;
        $options = 0, next  if $arg eq '--';
        push(@lhs, $arg), next  if $options;
        push(@rhs, $arg);
    }
    #warn "lhs/rhs: [@lhs] [@rhs]\n";

    push @errors, $self->parse_options(\@lhs);

    #fail unprocessed options; leave implicit non-options
    @$argv = ();
    foreach my $arg (@lhs) {
        if ($arg =~ /^--?\S/) {
            push @errors, "unknown or bad option '$arg'";
        }
        push @$argv, $arg;
    }
    push @$argv, @rhs;  #replace explicit non-options

    #save remaining ARGV for posterity
    $self->{'param'}->{'argvfiles'} = [@$argv];

    #instantiate global parameter object
    unless (@errors) {
        new Bio::MView::Option::Parameters($self->{'param'});
    }

    return @errors;
}

sub get_parameters { return $_[0]->{'param'} }

sub usage {
    my $self = shift;
    my $s = $Bio::MView::Option::Options::Header;  #already loaded
    my $prog = $self->{'param'}->{'prog'};
    $s =~ s/<PROG>/$prog/;

    foreach my $grp (@{$self->{'groups'}}) {
        next  unless exists $grp->{'header'};  #hidden group

        $s .= $grp->{'header'} . "\n";

        foreach my $opt (@{$grp->{'options'}}) {
            $s .= $self->get_option_usage($opt);
        }

        $s .= "\n";
    }

    return $s;
}

######################################################################
# private class methods
######################################################################
sub str {
    return '<UNDEF>'  unless defined $_[0];
    return "[@{$_[0]}]"  if ref $_[0] eq 'ARRAY';
    return $_[0];
}

sub load_library {
    my $library = $_[0];
    $library =~ s/::/\//g;
    require "$library.pm";
}

sub load_scalar {
    my ($lib, $var) = @_;
    return eval '$' . $lib . '::' . $var;
}

######################################################################
# private methods
######################################################################
sub load_types {
    my ($self, $lib) = @_;
    my $var = load_scalar($lib, 'Types');
    foreach my $type (@$var) {
        my $tn = $type->{'type'};
        if (exists $self->{'types'}->{$tn}) {
            die "Options: type '$tn' already exists\n";
        }
        $self->{'types'}->{$tn} = $type;
    }
    return $self;
}

sub load_option_groups {
    my ($self, $lib, $seen) = @_;
    my $var = load_scalar($lib, 'Options');
    foreach my $grp (@$var) {
        my $gn = $grp->{'group'};
        if (exists $seen->{$gn}) {
            die "Options: option group '$gn' already exists\n";
        }
        if (! exists $grp->{'options'}) {
            die "Options: option group '$gn' contains no options\n";
        }
        push @{$self->{'groups'}}, $grp;
        $seen->{$gn} = $grp;
        $self->load_options($grp, $gn);
    }
    return $self;
}

sub load_options {
    my ($self, $grp, $gn) = @_;
    foreach my $opt (@{$grp->{'options'}}) {
        my $on = $gn . '::' . $opt->{'option'};
        if (exists $self->{'options'}->{$on}) {
            die "Options: option '$on' already exists\n";
        }
        $opt->{'group'} = $gn;  #insert extended option name
        $self->{'options'}->{$on} = $opt;
        $self->initialise_option($opt);
    }
    return $self;
}

sub initialise_option {
    my ($self, $opt) = @_;

    my $gn = $opt->{'group'};
    my $on = $opt->{'option'};

    if (! exists $opt->{'type'}) {
        die "Options: option '$gn.$on' has no type\n";
    }

    if (! $opt->{'type'}) {
        die "Options: option '$gn.$on' has empty type\n";
    }

    my $type = $self->get_type_record($opt);

    if (! defined $type) {
        die "Options: option '$gn.$on' has unbound type\n";
    }

    #set an option ordering
    if (! exists $opt->{'order'}) {
        $opt->{'order'} = keys %{$self->{'options'}};
    }

    #set a default value, unless set already
    if (! exists $opt->{'default'}) {

        if (! exists $type->{'default'}) {
            warn "Options: no default for option '$gn.$on'\n";
            die  "Options: no default for type '$type->{'type'}'\n";
        }

        $opt->{'default'} = $type->{'default'};
    }

    #set the default action, unless one specified
    $opt->{'action'} = undef  unless exists $opt->{'action'};

    #set the option value from the default: will be overwritten from CL
    $opt->{'option_value'} = $opt->{'default'};

    #set the associated parameter name, unless set already
    $opt->{'param'} = $opt->{'option'}  unless exists $opt->{'param'};

    #set the default parameter value
    my @errors = $self->test_and_set_option($opt);

    if (@errors) {
        warn "Options: option '$gn.$on' initialisation errors:\n";
        foreach my $e (@errors) {
            warn "Options: $e\n";
        }
    }

    #dump_option($opt);

    return $self;
}

sub update_option {
    my ($self, $on, $ov) = @_;

    return  unless exists $self->{'options'}->{$on};

    my $opt = $self->{'options'}->{$on};

    $opt->{'option_value'} = $ov;

    $self->test_and_set_option($opt);
}

sub get_type_record {
    my ($self, $key) = @_;
    while (ref $key) {  #search for a base type
        $key = $key->{'type'};
    }
    return $self->{'types'}->{$key}  if exists $self->{'types'}->{$key};
    return undef;
}

sub get_attr_string {
    my ($self, $opt, $attr) = @_;
    return str($opt->{$attr})   if exists $opt->{$attr};
    my $type = $self->get_type_record($opt);
    return str($type->{$attr})  if exists $type->{$attr};
    return '';
}

sub get_label_string {
    my ($self, $opt) = @_;
    return $self->get_attr_string($opt, 'label');
}

sub get_usage_string {
    my ($self, $opt) = @_;
    return $self->get_attr_string($opt, 'usage');
}

sub get_values_string {
    my ($self, $opt) = @_;
    return $self->get_attr_string($opt, 'values');
}

sub get_default_string {
    my ($self, $opt) = @_;
    return $self->get_attr_string($opt, 'default');
}

sub get_option_usage {
    my ($self, $opt) = @_;
    my $name    = $opt->{'option'};
    my $type    = $self->get_label_string($opt);
    my $usage   = $self->get_usage_string($opt);
    my $values  = $self->get_values_string($opt);
    my $default = $self->get_default_string($opt);
    my $s = sprintf("  -%-20s %s", "$name $type", $usage);
    $s .= " {$values}"   if $values  ne '';
    $s .= ".";
    $s .= " [$default]"  if $default ne '';
    $s .= "\n";
    return $s;
}

sub get_option_format {
    my ($self, $opt) = @_;
    my $on = $opt->{'option'};

    #flag: no argument
    return $on  if $opt->{'type'} =~/^flag/;

    #option: takes argument
    return "$on=s";
}

sub test_type {
    my ($self, $opt, $on, $pv, $e) = @_;

    return  unless defined $pv;

    my $type = $self->get_type_record($opt);
    my $test = $type->{'test'};

    if (ref $test eq 'CODE') {
        $pv = &{$test}($self, $on, $pv, $e);
    }

    return $pv;
}

sub perform_action {
    my ($self, $opt, $pn, $pv, $e) = @_;

    #only look in $opt; search $typet currently not needed
    my $action = $opt->{'action'};

    return  unless ref $action eq 'CODE';
    return  unless defined $pv;

    &{$action}($self, $pn, $pv, $e);
}

sub test_and_set_option {
    my ($self, $opt) = @_;
    my @errors = ();

    #update the parameter hash with this option
    my $on = $opt->{'option'};
    my $ov = $opt->{'option_value'};

    #type tests and simple parameter conversion
    my $pn = $opt->{'param'};
    my $pv = $self->test_type($opt, $on, $ov, \@errors);

    #update the parameter hash
    $self->{'param'}->{$pn} = $pv;

    #immediately perform associated action if any
    if (defined $opt->{'action'}) {
        $self->perform_action($opt, $pn, $pv, \@errors);
    }

    if ($DEBUG) {
        printf STDERR "opt:%15s => %-10s    par:%15s => %-10s\n",
            $opt->{'option'}, str($opt->{'option_value'}),
            $opt->{'param'},  str($self->{'param'}->{$opt->{'param'}});
    }
    return @errors;
}

sub parse_options {
    my ($self, $argv) = @_;

    my @errors = (); my $scan = {};

    #process options in specific order; don't abort on error:
    foreach my $opt ($self->option_order) {

        my $format = $self->get_option_format($opt);

        next  unless $format;

        GetOptionsFromArray($argv, $scan, $format);

        #map { print STDERR "$_ => $scan->{$_}\n" } %$scan;

        my $on = $opt->{'option'};

        next  unless exists $scan->{$on};

        #save scanned value
        $opt->{'option_value'} = $scan->{$on};

        push @errors, $self->test_and_set_option($opt);
    }

    return @errors;
}

sub option_order {
    my $self = shift;
    return sort { $a->{'order'} <=> $b->{'order'} }
        values %{$self->{'options'}};
}

######################################################################
# debug
######################################################################
sub dump_item {
    my ($item, $stm) = (shift, shift);
    foreach my $key (@_) {
        next  unless exists $item->{$key};
        my $val = $item->{$key};
        $val = '<UNDEF>'  unless defined $val;
        $val = '<CODE>'     if ref $val eq 'CODE';
        $val = "''"       if length $val < 1;
        print $stm sprintf "%20s => %s\n", $key, $val;
    }
    print $stm "\n";
}

sub dump_type {
    my ($type, $stm) = (@_, \*STDERR);
    my @fields = qw(type label usage default test);
    dump_item($type, $stm, @fields);
}

sub dump_option {
    my ($opt, $stm) = (@_, \*STDERR);
    my @fields = qw(group option option_value usage default type test
                    label param);
    dump_item($opt, $stm, @fields);
}


###########################################################################
1;
