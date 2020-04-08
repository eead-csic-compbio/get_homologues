# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Util::System;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(stacktrace vmstat is_unix is_dos);

sub stacktrace {
    warn "Stack Trace:\n"; my $i = 0;
    my @calls = caller($i++);
    my ($file, $line, $func) = ($calls[1], $calls[2], $calls[3]);
    while ( @calls = caller($i++) ){
        #func is one ahead
        warn $file . ":" . $line . " in function " . $calls[3] . "\n";
        ($file, $line, $func) = ($calls[1], $calls[2], $calls[3]);
    }
}

#Linux only
sub vmstat {
    my ($s) = (@_, '');
    local ($_, *TMP);
    if (open(TMP, "cat /proc/$$/stat|")) {
        $_ = <TMP>; my @ps = split /\s+/; close TMP;
        print sprintf "VM: %8luk  $s\n", $ps[22] >> 10;
    } else {
        print sprintf "VM: -  $s\n";
    }
}

sub is_unix {
    return 1
        if grep {/$^O/i} qw(aix android bsdos beos bitrig dgux dynixptx cygwin
                            darwin dragonfly freebsd gnukfreebsd haiku hpux
                            interix irix linux machten midnightbsd minix
                            mirbsd netbsd next nto openbsd qnx sco solaris);
    return 0;
}

sub is_dos {
    return 1
        if grep {/$^O/i} qw(dos MSWin32);
    return 0;
}

###########################################################################
1;
