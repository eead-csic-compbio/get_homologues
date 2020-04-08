# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::Util::Object;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(dump_self dump_hash);

my $PRETTY_PRINT_STREAM = \*STDOUT;
my $MIN_KEY_WIDTH = 16;

###########################################################################
# public methods
###########################################################################
# warn with error string
sub warn {
    my $self = shift;
    my $s = $self->make_message_string('Warning', @_);
    warn "$s\n";
}

# exit with error string
sub die {
    my $self = shift;
    my $s = $self->make_message_string('Died', @_);
    die "$s\n";
}

# return basename component of class
sub basename_class {
    my @id = split("::", ref($_[0]));
    return pop @id;
}

# pretty-print all instance variables on object, or supplied list
sub examine {
    print $PRETTY_PRINT_STREAM dump_self(@_);
    return $_[0];
}

###########################################################################
# exported statics
###########################################################################
# return pretty-printed object contents by given ordered keys, or all sorted
sub dump_self { return "Class $_[0]\n" . _dump_body(@_) }

# return pretty-printed hash contents by given ordered keys, or all sorted
sub dump_hash { return "$_[0]\n" . _dump_body(@_) }

###########################################################################
# protected functions
###########################################################################
# generate common error message string
sub make_message_string {
    my ($self, $prefix) = (shift, shift);
    my $s = "$prefix ";
    $s .= ref($self) ? ref($self) : $self;
    $s .= ": " . _args_as_string(@_)  if @_;
    return $s;
}

###########################################################################
# private functions
###########################################################################
use Bio::Util::Regexp;
use Bio::Util::Math qw(max);

# return space-separated string of concatenate arguments with undefined values
# interpolated as 'undef'
sub _args_as_string {
    my @tmp = ();
    foreach my $a (@_) {
        push @tmp, (defined $a ? $a : 'undef');
    }
    my $s = join(" ", @tmp);
    chomp $s;
    return $s;
}

# return pretty-printed hash contents by given ordered keys, or all sorted
sub _dump_body {
    my $hash = shift;

    sub maxlen {
        my $w = $MIN_KEY_WIDTH;
        foreach my $key (@_) {
            $w = max($w, length(defined $key ? $key : '<NOEXIST>'));
        }
        return $w;
    }

    sub layout {
        my $w = shift; return sprintf("%${w}s => %s\n", @_);
    }

    push @_, sort keys %$hash  unless @_;

    my $w = maxlen(@_); return ''  unless $w > 0;
    my $s = '';

    foreach my $key (@_) {
        if (exists $hash->{$key}) {
            my $val = $hash->{$key};
            if (! defined $val) {
                $s .= layout($w, $key, 'undef');
                next;
            }
            my $ref = ref $val;
            if ($ref) {
                if ($ref eq 'ARRAY') {
                    $s .= layout($w, $key, "@[" . join(',', @$val)  . "]");
                    next;
                }
                if ($ref eq 'HASH') {
                    my @tmp = map { "$_:$val->{$_}" } sort keys %$val;
                    $s .= layout($w, $key, "%{" . join(',', @tmp) . "}");
                    next;
                }
                $s .= layout($w, $key, "<$ref>");  #other ref
                next;
            }
            # numeric
            if ($val =~ /^$RX_Sreal$/) {
                $s .= layout($w, $key, (length($val) > 0 ? $val : "'$val'"));
                next;
            }
            #string
            $s .= layout($w, $key, "'$val'");
        } else {
            $s .= layout($w, $key, 'NOEXIST');
        }
    }
    return $s;
}

######################################################################
1;
