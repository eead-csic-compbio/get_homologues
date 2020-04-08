# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA1::fasta;

use Bio::Parse::Format::FASTA1;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA1);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA1::fasta::HEADER;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::HEADER);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'query'} = '';

    while (defined ($line = $scan->read_line)) {

        if ($line =~ /\s*(version\s+(\S+).*)/) {
            $self->{'full_version'} = $1;
            $self->{'version'}      = $2;
        }

        if ($line =~ /(\S+)\s*,\s+(\d+)\s+(?:aa|nt)/o) {

            $self->test_args(\$line, $1, $2);

            (
             $self->{'queryfile'},
             $self->{'length'},
            ) = ($1, $2);

            next;
        }

        if ($line =~ /^\s*>(\S+)\s*\:\s*\d+\s+(?:aa|nt)/) {
            $self->test_args(\$line, $1);
            $self->{'query'} = clean_identifier($1);
            next;
        }

        if ($line =~ /^(\d+)\s+residues\s+in\s+(\d+)\s+sequences/) {

            $self->test_args(\$line, $1,$2);

            (
             $self->{'residues'},
             $self->{'sequences'},
            ) = ($1, $2);

            next;
        }

        #ignore any other text

    }

    if (! defined $self->{'full_version'} ) {
        #can't determine version: hardwire one!
        $self->{'full_version'} = 'looks like FASTA 1';
        $self->{'version'}      = '1?';
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA1::fasta::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #ranked search hits
    while (defined ($line = $scan->read_line)) {

        next    if $line =~ /$Bio::Parse::Format::FASTA1::RANK_START/o;

        if ($line =~ /^
            \s*
            (\S+)                #id
            \s+
            (.*)                 #description
            \s+
            (\d+)                #initn
            \s+
            (\d+)                #init1
            \s+
            (\d+)                #opt
            \s*
            $/xo) {

            $self->test_args(\$line, $1,$2,$3,$4,$5);

            push @{$self->{'hit'}},
            {
             'id'    => clean_identifier($1),
             'desc'  => $2,
             'initn' => $3,
             'init1' => $4,
             'opt'   => $5,
            };

            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA1::NULL/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA1::fasta::TRAILER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA1::fasta::MATCH;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA1::fasta::MATCH::SUM;

use Bio::Parse::Strings qw(strip_english_newlines clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::MATCH::SUM);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $line = $scan->read_line;

    if ($line =~ /^(\S+)\s+(.*)\s+(\d+)\s+(\d+)\s+(\d+)\s*$/) {

        $self->test_args(\$line, $1, $2, $3, $4, $5);    #ignore $2

        (
         $self->{'id'},
         $self->{'desc'},
         $self->{'initn'},
         $self->{'init1'},
         $self->{'opt'},
        ) = (clean_identifier($1), strip_english_newlines($2), $3, $4, $5);

    }

    $line = $scan->read_line;

    if ($line =~ /^\s*($RX_Ureal)\% identity in (\d+) (?:aa|nt) overlap\s*$/) {
        $self->test_args(\$line, $1, $2);
        (
         $self->{'id_percent'},
         $self->{'overlap'},
        )= ($1,$2);
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA1::fasta::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA1::MATCH::ALN);


###########################################################################
1;
