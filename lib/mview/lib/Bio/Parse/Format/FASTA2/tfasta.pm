# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA2::tfasta;

use Bio::Parse::Format::FASTA2;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA2);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA2::tfasta::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA2::HEADER);


###########################################################################
package Bio::Parse::Format::FASTA2::tfasta::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::RANK);

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
           \((\S)\)                #frame
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

            $self->test_args(\$line, $1,$2,$3,$4,$5,$6,$7,$8);

            push(@{$self->{'hit'}},
                 {
                  'id'     => clean_identifier($1),
                  'desc'   => $2,
                  'frame'  => Bio::Parse::Format::FASTA::parse_frame($3),
                  'orient' => Bio::Parse::Format::FASTA::parse_orient($3),
                  'initn'  => $4,
                  'init1'  => $5,
                  'opt'    => $6,
                  'zscore' => $7,
                  'expect' => $8,
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
package Bio::Parse::Format::FASTA2::tfasta::TRAILER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA2::tfasta::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA2::tfasta::MATCH::SUM;

use Bio::Parse::Strings qw(strip_english_newlines clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH::SUM);

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
        (?:Frame\:\s*(\S+))?   #frame
        \s*
        init1\:\s*(\S+)        #init1  REVERSE of fasta ordering WRONG!
        \s*
        initn\:\s*(\S+)        #initn  REVERSE of fasta ordering WRONG!
        \s*
        opt\:\s*(\S+)          #opt
        \s*
        z-score\:\s*(\S+)      #z
        \s*
        E\(\)\:\s*(\S+)        #E
        \s*
        $/xo) {

        $self->test_args(\$line,$2,$3,$4,$5,$6);

        (
         $self->{'frame'},
         $self->{'orient'},
         $self->{'initn'},     #correct REVERSAL
         $self->{'init1'},     #correct REVERSAL
         $self->{'opt'},
         $self->{'zscore'},
         $self->{'expect'},
        ) = (
            Bio::Parse::Format::FASTA::parse_frame($1),
            Bio::Parse::Format::FASTA::parse_orient($1),
            $2, $3, $4, $5, $6,
        );
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
package Bio::Parse::Format::FASTA2::tfasta::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA2::MATCH::ALN);


###########################################################################
1;
