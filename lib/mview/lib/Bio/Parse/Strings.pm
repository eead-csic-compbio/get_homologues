# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Parse::Strings;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(
    strip_english_newlines
    strip_leading_space
    strip_trailing_space
    strip_leading_identifier_chars
    strip_trailing_identifier_chars
    clean_identifier
);

# return $s less newlines while attempting to assemble hyphenated words and
# excess white space split over multiple lines correctly.
sub strip_english_newlines {
    my $s = shift;

    #multiple hyphens look like 'SO4--' or long dashes - word break with space
    $s =~ s/(--)\s*\n+\s*/$1 /sg;

    #single hyphens look like hyphenated words - join at the hyphen
    $s =~ s/(-)\s*\n+\s*/$1/sg;

    #remaining newlines - word break with space
    $s =~ s/\s*\n+\s*/ /sg;

    #skip trailing white space added after last newline was removed
    $s =~ s/\s+$//s;

    return $s;
}

sub strip_leading_space {
    my $s = shift;
    $s =~ s/^[ \t]+//;
    return $s;
}

sub strip_trailing_space {
    my $s = shift;
    $s =~ s/[ \t]+$//;
    return $s;
}

sub strip_leading_identifier_chars {
    my $s = shift;
    $s =~ s/^(\/:|>)//;  #leading "/:" or ">"
    return $s;
}

sub strip_trailing_identifier_chars {
    my $s = shift;
    $s =~ s/[;:|,.-]+$//;  #default trailing PDB chain is often '_'
    return $s;
}

sub clean_identifier {
    my $s = shift;
    $s = strip_leading_identifier_chars($s);
    $s = strip_trailing_identifier_chars($s);
    return $s;
}

###########################################################################
1;
