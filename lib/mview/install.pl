#!/usr/bin/env perl

# Copyright (C) 2019-2020 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
# application to be installed (UNIX syntax; will be modified for DOS)
my $TARGET_NAME = "MView";
my $TARGET_ENV  = "MVIEW_HOME";
my $TARGET_EXE  = "bin/mview";
my $DRIVER      = "mview";
###########################################################################

require 5.004;  # ancient: same as MView

use strict;
use warnings;

my $VERSION  = "1.1"; # installer version
my $TESTFILE = "";    # tmpfile, unlink after Ctrl-C

sub abort {
    unlink($DRIVER)    if -f $DRIVER;
    unlink($TESTFILE)  if -f $TESTFILE;
    die "\nExited.\n";
}

# trap 1 2 3 15
$SIG{'HUP'}  = \&abort;
$SIG{'INT'}  = \&abort;
$SIG{'QUIT'} = \&abort;
$SIG{'TERM'} = \&abort;

sub show_admin_warning {
    print STDERR <<EOT;
###########################################################################

                     $TARGET_NAME installer

  *********************************************************
  **                                                     **
  **  WARNING: You have administrator/root permissions!  **
  **                                                     **
  *********************************************************

EOT
}

sub show_path_warning {
    print STDERR <<EOT;
  *********************************************************
  ** Note: That folder does not exist yet.               **
  **                                                     **
  ** Don't forget to add it manually to PATH.            **
  **                                                     **
  **                                                     **
  ** For another choice press Ctrl-C, then start again.  **
  *********************************************************

EOT
}#'

sub show_preamble {
    print STDERR <<EOT;
###########################################################################

                         $TARGET_NAME installer

The installation requires a folder to contain the driver script.

That folder must be (1) writable by the installer and (2) on your PATH if
installing a personal copy, or on the shared PATH for site installations.

You can accept the suggested default, choose a folder from the following
PATH list, or specify another folder that you will add to PATH later:

EOT
}

sub show_summary {
    my $installer = shift;
    my $installdir  = $installer->installdir();
    my $bindir      = $installer->bindir();
    my $driver_path = $installer->driver_path();
    my $short_name  = $installer->driver_short_name();
    my $long_name   = $installer->driver_long_name();

    print STDERR <<EOT;

###########################################################################

Installation details:

Installed under:  $installdir
PATH entry:       $bindir           (add to PATH if missing)
Program:          $driver_path

Installation complete.

Try running the program with:

  $short_name -help

If PATH is not set properly this should work:

  $long_name -help

###########################################################################
EOT
}

sub info    { print STDERR "@_\n"; }
sub warning { print STDERR " + @_\n"; }
sub error   { print STDERR "\nERROR: @_\n"; }
sub prompt  { print STDERR "\n>> @_ "; }

sub pause {
    prompt "Press the 'return' or 'enter' key to continue (Ctrl-C to exit).";
    <STDIN>;
    info;
}

sub pause_before_exit {
    prompt "Press the 'return' or 'enter' key to finish.";
    <STDIN>;
}

###########################################################################
package Installer;

use Cwd qw(getcwd);
use File::Basename;
use File::Copy qw(copy);
use File::Path qw(mkpath);
use File::Spec;
use POSIX qw(strftime);

sub new {
    my ($cls) = @_;
    my $self = {};

    $self->{'ADMIN'}      = undef;
    $self->{'INSTALLDIR'} = getcwd();
    $self->{'BINDIR'}     = undef;
    $self->{'SAVEPATH'}   = [];

    bless $self, $cls;
}

sub get_system {
    if (grep {/$^O/i} qw(dos MSWin32)) {
        return "DOS";
    }
    if (grep {/$^O/i} qw(aix android bsdos beos bitrig dgux
        dynixptx cygwin darwin dragonfly freebsd gnukfreebsd haiku hpux interix
        irix linux machten midnightbsd minix mirbsd netbsd next nto openbsd qnx
        sco solaris)) {
        return "UNIX";
    }
    return $^O;
}

sub _get_timestamp {
    my $self = shift;
    # Mon Sep 30 17:31:18 CEST 2019
    strftime "%a %b %d %H:%M:%S %Z %Y", localtime;
}

# _make_dir dirpath
sub _make_dir {
    my ($self, $dir) = @_;
    ::warning "Making directory '$dir'";
    return mkpath($dir);
}

# _make_symlink target linkname
sub _make_symlink () {
    my ($self, $target, $link) = @_;
    ::warning "Making symlink from '$link' to '$target'";
    return symlink($target, $link);
}

# _make_copy source destination
sub _make_copy {
    my ($self, $source, $destination) = @_;
    ::warning "Making copy from '$source' to '$destination'";
    return copy($source, $destination);
}

# _remove_file filepath
sub _remove_file {
    my ($self, $file) = @_;
    ::warning "Deleting '$file'";
    if (!( -f $file or -l $file )) {
        ::warning "Name '$file' is not a file or symlink";
        return 1;
    }
    return unlink($file);
}

sub make_bindir {
    my $self = shift;
    my $bindir = $self->{'BINDIR'};
    if ( ! -e $bindir ) {
        if ($self->_make_dir($bindir) < 1) {
            ::warning "Can't make directory '$bindir'";
        } else {
            ::warning "Created directory '$bindir'";
        }
    }
}

sub list_dirs {
    my $self = shift;
    if (@{$self->{'SAVEPATH'}}) {
        foreach my $p (@{$self->{'SAVEPATH'}}) {
            ::info "  $p";
        }
    } else {
        ::info "No available paths";
    }
}

sub cleanup {
    my $self = shift;
    unlink($self->{'DRIVER'})  if -f $self->{'DRIVER'};
}

sub in_package_dir {
    return 0  unless -f basename($0);
    return 0  unless -f $_[0]->{'TARGET_EXE'};
    return 1;
}

sub installdir { return $_[0]->{'INSTALLDIR'} }

sub bindir { return $_[0]->{'BINDIR'} }

sub driver_path {
    return File::Spec->catfile($_[0]->{'BINDIR'}, $_[0]->{'DRIVER'});
}

sub driver_short_name { return $_[0]->{'DRIVER'} }

sub driver_long_name {
    File::Spec->catfile($_[0]->{'BINDIR'}, $_[0]->driver_short_name());
}

sub set_bindir {
    my ($self, $dir) = @_;
    $self->{'BINDIR'} = $self->_expand_path($dir);
}

sub guess_bindir {
    my $self = shift;

    return $self->{'BINDIR'}  if defined $self->{'BINDIR'};  # already set

    $self->{'SAVEPATH'} = [ $self->_get_writable_paths() ];

    return $self->_guess_admin_bindir(@{$self->{'SAVEPATH'}})  if $self->is_admin();
    return $self->_guess_user_bindir(@{$self->{'SAVEPATH'}});
}

sub choose_bindir {
    my $self = shift;
    my ($bindir, $loc) = ($self->{'BINDIR'}, "");

    while (1) {
        ::prompt "Choose a directory for the program [$bindir]";
        my $raw = <STDIN>; chomp $raw;
        $loc = $self->_expand_path($raw);
        ::info;

        if ( $loc eq "" and $bindir ne "") {
            # accept default BINDIR
            $loc = $bindir;
            last;
        }
        if ( $loc eq "" ) {
            ::warning "Directory must have a name";
            next;
        }
        if ( $loc eq $self->{'DRIVER'} ) {
            ::warning "That name is reserved - please use another one";
            next;
        }
        if ( -e $loc and ! -d $loc ) {
            ::warning "Name '$raw' is not a directory";
            next;
        }
        if ( ! $self->_path_writable($loc) ) {
            ::warning "Directory path '$raw' is not writable";
            next;
        }
        # choice accepted - dir already exists or will be created by user
        last;
    }

    $bindir = $loc;

    ::show_path_warning()  unless -e $bindir;

    return $self->{'BINDIR'} = $bindir;
}

# return 1 if the path is potentially writable, 0 otherwise
sub _path_writable {
    my ($self, $path) = @_;

    while ($path) {
        if (-e $path and ! -d $path) {
            return 0;  # exists but not a directory
        }
        if (-e $path and $self->_is_writable_dir($path)) {
            return 1;  # exists and can write into it
        }
        if (-e $path and ! $self->_is_writable_dir($path)) {
            return 0;  # can't write below path element
        }
        # path element doesn't exist yet - potentially writable;
        # go up one level
        my ($v,$d,$f) = File::Spec->splitpath($path);

        return 1  if $d eq "";  # top-level of relative path

        $path = File::Spec->join($v, $d);
    }
    return 1;
}

sub install_driver {
    my $self = shift;

    return  if $self->{'BINDIR'} eq ".";

    my ($destination, $file) = ($self->{'BINDIR'}, $self->{'DRIVER'});

    my $source = "$file";
    my $target = File::Spec->catfile($destination, $file);

    if ( -e $target ) {
        ::warning "Replacing existing '$target'";
        if ($self->_remove_file($target) < 1) {
            $self->error("Attempt to delete old '$target' failed - exiting!");
        }
    }

    if ($self->_make_copy($source, $target) < 1) {
        $self->error("Attempt to copy '$source' into '$target' failed - exiting!");
    }

    $self->subclass_install_driver($target);
}

sub error {
    shift;
    ::error @_;
    exit 1;
}

sub subclass_install_driver {}

###########################################################################
package UnixInstaller;

our @ISA = qw(Installer);

sub new {
    my $cls = shift;
    my $self = new Installer(@_);

    $self->{'TARGET_EXE'} = $TARGET_EXE;
    $self->{'DRIVER'}     = $DRIVER;

    umask 022;

    bless $self, $cls;
}

sub is_dos  { 0 }

sub is_admin {
    my $self = shift;

    return $self->{'ADMIN'}  if defined $self->{'ADMIN'};

    my $user;
    if ($ENV{'USER'} eq "root") {
        $user = 0;
    }
    elsif ( -e "/usr/bin/id" ) {
        $user = `/usr/bin/id -u`;
    }
    elsif ( -e "/bin/id" ) {
        $user = `/bin/id -u`;
    }
    else {
        $user = `id -u`;
    }

    return $self->{'ADMIN'} = ($user == 0);
}

sub _is_writable_dir {
    my ($self, $dir) = @_;
    return 1  if -w $dir;
    return 0;
}

sub _expand_path {
    my ($self, $path) = @_;
    $path =~ s{ ^~([^/]*) } { $self->_get_home_dir($1) }ex;
    return $path;
}

sub _get_home_dir {
    my ($self, $s) = @_;
    return (getpwnam($s))[7]  if @_ and $s;
    return $ENV{HOME}         if exists $ENV{HOME};
    return $ENV{LOGDIR}       if exists $ENV{LOGDIR};
    return (getpwuid($<))[7];
}

sub _get_writable_paths {
    my $self = shift;

    my $paths = $ENV{'PATH'};
    my (%seen, @list) = ();

    foreach my $p (split(':', $paths)) {
        $p = $self->_expand_path($p);
        next  unless -w $p;  # writable
        push @list, $p  unless exists $seen{$p};
        $seen{$p} = 1;
    }

    return @list;
}

sub _guess_admin_bindir {
    my $self = shift;

    # prefer common paths in this order
    foreach my $dir (
        '/usr/local/bin',
        '/opt/bin',
    ) {
        if (my @tmp = grep { /^$dir$/ } @_) {
            return $self->{'BINDIR'} = $tmp[0];
      }
    }

    return $self->{'BINDIR'} = "/usr/local/bin";  # force default
}

sub _guess_user_bindir {
    my $self = shift;

    my $HOME = $self->_get_home_dir();

    # prefer paths in this order
    foreach my $dir (
        "$HOME/bin",
    ) {
        if (my @tmp = grep { /^$dir$/ } @_) {
            return $self->{'BINDIR'} = $tmp[0];
        }
    }

    return $self->{'BINDIR'} = "$HOME/bin";  # force default
}

sub make_driver {
    my $self = shift;

    my $file = $self->{'DRIVER'};
    my $date = $self->_get_timestamp();
    my $installdir = $self->{'INSTALLDIR'};
    my $executable = $self->{'TARGET_EXE'};

    open(my $fh, ">", $file);
    print $fh <<EOT;
#!/bin/sh
# $TARGET_NAME driver
# Version: $VERSION
# Generated: $date

$TARGET_ENV=$installdir; export $TARGET_ENV
PROGRAM=\$$TARGET_ENV/$executable
PROG=\`basename \$0\`

# echo $TARGET_ENV=\$$TARGET_ENV  1>&2
# echo PROGRAM=\$PROGRAM  1>&2

if [ ! -f \$PROGRAM ]; then
    echo "\$PROG: Can't find program '\$PROGRAM'"
    exit 1
fi
if [ ! -x \$PROGRAM ]; then
    echo "\$PROG: Program '\$PROGRAM' is not executable"
    exit 1
fi

exec \$PROGRAM "\$@"
EOT
    close $fh;
}

# override
sub subclass_install_driver {
    my ($self, $target) = @_;
    my $mode = 0755;
    if (chmod($mode, $target) < 1) {
        $self->error("Attempt to set execute permissions on '$target' failed - exiting!");
    }
}

###########################################################################
package DosInstaller;

our @ISA = qw(Installer);

sub new {
    my $cls = shift;
    my $self = new Installer(@_);

    #change from defaults
    $self->{'INSTALLDIR'} = join('\\', split('/', $self->{'INSTALLDIR'}));
    $self->{'TARGET_EXE'} = join('\\', split('/', $TARGET_EXE));
    $self->{'DRIVER'}     = "$DRIVER.bat";

    # set global for trap
    $DRIVER = $self->{'DRIVER'};

    bless $self, $cls;
}

sub is_dos { 1 }

sub driver_short_name {
    my $self = shift;
    my $tmp = $self->{'DRIVER'};
    $tmp =~ s/\.bat$//i;  # strip batch file extension
    return $tmp;
}

sub is_admin {
    my $self = shift;

    return $self->{'ADMIN'}  if defined $self->{'ADMIN'};

    $self->{'ADMIN'} = 0;
    $self->{'ADMIN'} = 1  if system("NET SESSION >NUL 2>&1") == 0;

    return $self->{'ADMIN'};
}

# ugly hack
sub _is_writable_dir {
    my ($self, $dir) = @_;

    # may seem 'writable' yet not writable
    return 0  unless -w $dir;

    # try to write a temporary file

    # set global for trap
    $TESTFILE = join('\\', $dir, "mview_$$");
    open(my $fh, ">", $TESTFILE) or return 0;
    close $fh;
    unlink $TESTFILE;

    return 1;
}

sub _expand_path {
    my ($self, $path) = @_;
    $path =~ s{ ^%UserProfile% } { $self->_get_home_dir() }iex;
    $path =~ s{ ^~ }             { $self->_get_home_dir() }iex;  #unix tilde
    return $path;
}

sub _get_home_dir {
    my $self = shift;
    return $ENV{'UserProfile'}  if exists $ENV{'UserProfile'};
    return "";
}

sub _get_writable_paths {
    my $self = shift;

    my $paths = $ENV{'PATH'};
    my (%seen, @list) = ();

    foreach my $p (split(';', $paths)) {
        next  unless $self->_is_writable_dir($p);
        push @list, $p  unless exists $seen{uc $p}; # case insensitive
        $seen{uc $p} = 1;
    }

    return @list;
}

sub _guess_admin_bindir {
    my $self = shift;

    # prefer common paths in this order
    foreach my $path (
        'C:\strawberry\perl\site\bin',
        'C:\perl\site\bin',
        'C:\perl\bin',
        'C:\bin',
    ) {
        my $dir = quotemeta($path);

        if (my @tmp = grep { /^$dir$/i } @_) {
            return $self->{'BINDIR'} = $tmp[0];
        }
    }

    return $self->{'BINDIR'} = 'C:\bin';  # force default
}

sub _guess_user_bindir {
    my $self = shift;

    my $HOME = $self->_get_home_dir();

    # prefer paths in this order
    foreach my $path (
        "$HOME\\bin",
    ) {
        my $dir = quotemeta($path);

        if (my @tmp = grep { /^$dir$/i } @_) {
            return $self->{'BINDIR'} = $tmp[0];
        }
    }

    return $self->{'BINDIR'} = "$HOME\\bin"  # force default
}

sub make_driver {
    my $self = shift;

    my $file = $self->{'DRIVER'};
    my $date = $self->_get_timestamp();
    my $installdir = $self->{'INSTALLDIR'};
    my $executable = $self->{'TARGET_EXE'};

    open(my $fh, ">", $file);
    print $fh <<EOT;
\@echo off
rem $TARGET_NAME driver
rem Version: $VERSION
rem Generated: $date

set $TARGET_ENV=$installdir
set PROGRAM=%$TARGET_ENV%\\$executable
set PROG=%~nx0

rem echo "$TARGET_ENV=%$TARGET_ENV%"  1>&2
rem echo "PROGRAM=%PROGRAM%"  1>&2

if not exist "%PROGRAM%" (
    echo "%PROG%: Can't find program '%PROGRAM%'"
    exit /b 1
)

perl "%PROGRAM%" %*
EOT
    close $fh;
}

# override
sub error {
    shift;
    ::error @_;
    ::pause_before_exit();  # hold any terminal window open
    exit 1;
}

###########################################################################
package Installer;

sub installer {
    my $system = get_system();
    return new UnixInstaller()  if $system eq "UNIX";
    return new DosInstaller()   if $system eq "DOS";
    ::error "I do not recognise the operating system '$^O' - exiting!";
    ::pause_before_exit();  # hold any terminal window open
    exit 1;
}

###########################################################################
package main;

my $installer = Installer::installer();

if (! $installer->in_package_dir()) {
    $installer->error("You must change into the unpacked directory first - exiting!");
}

$installer->cleanup();  # any previous attempt

if ($installer->is_admin()) {
    show_admin_warning();
    pause();
}

if (@ARGV) {
    # destination directory on command line
    my $bindir = shift;
    $installer->set_bindir($bindir);
    info;
    info "Installing driver script into '$bindir'";
    info;
}
else {
    # interactive choice
    my $bindir;
    show_preamble();
    $bindir = $installer->guess_bindir();
    $installer->list_dirs();
    info;
    info "The suggested default is [$bindir]";
    info;
    info "At any time, you may press (Ctrl-C) to exit.";
    $bindir = $installer->choose_bindir();
    info "About to install driver script into '$bindir'";
    pause()  if @ARGV < 1;
}

$installer->make_bindir();
$installer->make_driver();
$installer->install_driver();
$installer->cleanup();

show_summary($installer);
pause_before_exit()  if $installer->is_dos();
exit 0;
