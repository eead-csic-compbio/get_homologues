# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA2::fasta;

use Bio::Parse::Format::FASTA2;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA2);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA2::fasta::HEADER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2::HEADER);


###########################################################################
package Bio::Parse::Format::FASTA2::fasta::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #ranked search hits
    while (defined ($line = $scan->read_line)) {

        next    if $line =~ /$Bio::Parse::Format::FASTA2::RANK_START/o;

        if($line =~ /^
           \s*
           ([^\s]+)                #id
           \s+
           (.*)                    #description
           \s+
           (\d+)                   #initn
           \s+
           (\d+)                   #init1
           \s+
           (\d+)                   #opt
           \s+
           (\S+)                   #z-score
           \s+
           (\S+)                   #E(58765)
           \s*
           $/xo) {

            $self->test_args(\$line, $1,$2,$3,$4,$5,$6,$7);

            push(@{$self->{'hit'}},
                 {
                  'id'     => clean_identifier($1),
                  'desc'   => $2,
                  'initn'  => $3,
                  'init1'  => $4,
                  'opt'    => $5,
                  'zscore' => $6,
                  'expect' => $7,
                 });
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA2::NULL/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA2::fasta::TRAILER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA2::fasta::MATCH;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA2::fasta::MATCH::SUM;

use Bio::Parse::Strings qw(strip_english_newlines clean_identifier);

use vars qw(@ISA);
use Bio::Util::Regexp;

@ISA   = qw(Bio::Parse::Format::FASTA::MATCH::SUM);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $line = $scan->read_line;

    if ($line =~ /^
        >*
        (\S+)                      #id
        \s+
        (.*)                       #description
        \s+
        \(\s*(\d+)\s*(?:aa|nt)\)   #length
        \s*
        $/xo) {

        $self->test_args(\$line, $1, $2, $3);

        (
         $self->{'id'},
         $self->{'desc'},
         $self->{'length'},
        ) = (clean_identifier($1), strip_english_newlines($2), $3);
    } else {
        $self->warn("unknown field: $line");
    }

    $line = $scan->read_line;

    if ($line =~ /^
        initn\:\s*(\S+)        #initn
        \s*
        init1\:\s*(\S+)        #init1
        \s*
        opt\:\s*(\S+)          #opt
        \s*
        z-score\:\s*(\S+)      #z
        \s*
        E\(\)\:\s*(\S+)        #E
        \s*
        $/xo) {

        $self->test_args(\$line,$1,$2,$3,$4,$5);

        (
         $self->{'initn'},
         $self->{'init1'},
         $self->{'opt'},
         $self->{'zscore'},
         $self->{'expect'},
        ) = ($1,$2,$3,$4,$5);
    } else {
        $self->warn("unknown field: $line");
    }

    $line = $scan->read_line;

    if ($line =~ /^
        (?:Smith-Waterman\s+score:\s*(\d+);)?    #sw score
        \s*($RX_Ureal)%                          #percent identity
        \s*identity\s+in\s+(\d+)                 #overlap length
        \s+(?:aa|nt)\s+overlap
        \s*
        $/xo) {

        $self->test_args(\$line,$2,$3);

        (
         $self->{'score'},
         $self->{'id_percent'},
         $self->{'overlap'},
        ) = (defined $1?$1:0,$2,$3);
    } else {
        $self->warn("unknown field: $line");
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA2::fasta::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA2::MATCH::ALN);


###########################################################################
1;
