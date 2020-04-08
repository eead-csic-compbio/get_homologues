# Copyright (C) 2011-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA3X::tfasta;

use Bio::Parse::Format::FASTA3X;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA3X::tfasta::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X::HEADER);


###########################################################################
package Bio::Parse::Format::FASTA3X::tfasta::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #ranked search hits
    while (defined ($line = $scan->read_line)) {

        next    if $line =~ /$Bio::Parse::Format::FASTA3X::RANK_START/o;

        if($line =~ /^
           \s*
           (\S+)                #id
           \s+
           (.*)                 #desc (may be empty)
           \s+
           \(\s*(\d+)\)         #aa
           \s+
           \[(\S)\]             #frame
           \s+
           (\d+)                #initn
           \s+
           (\d+)                #opt
           \s+
           (\S+)                #bits
           \s+
           (\S+)                #E(205044)
           \s*
           $/xo) {


            $self->test_args(\$line, $1, $3,$4, $5,$6,$7,$8); #not $2

            push(@{$self->{'hit'}},
                 {
                  'id'     => clean_identifier($1),
                  'desc'   => $2,
                  'length' => $3,
                  'frame'  => Bio::Parse::Format::FASTA::parse_frame($4),
                  'orient' => Bio::Parse::Format::FASTA::parse_orient($4),
                  'initn'  => $5,
                  'init1'  => 0,
                  'opt'    => $6,
                  'bits'   => $7,
                  'zscore' => 0,
                  'expect' => $8,
                 });
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA3X::NULL/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA3X::tfasta::TRAILER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA3X::tfasta::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA3X::tfasta::MATCH::SUM;

use Bio::Parse::Strings qw(strip_english_newlines clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH::SUM);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    my $record = '';

    while (defined ($line = $scan->read_line(1))) {

        next  if $line =~ /^\s*$/;

        if ($line =~ /^>+/) {
            $record = $line;       #process this later
            next;
        }

        if ($line =~ /^
            Frame\:\s*(\S+)        #frame [fr]
            \s+
            initn\:\s*(\S+)        #initn
            \s+
            init1\:\s*(\S+)        #init1
            \s+
            opt\:\s*(\S+)          #opt
            \s+
            Z-score\:\s*(\S+)      #z
            \s+
            bits\:\s*(\S+)         #bits
            \s+
            E\((?:.*)\)\:\s*(\S+)  #E
            \s*
            $/xo) {

            $self->test_args(\$line,$1,$2,$3,$4,$5,$6,$7);

            (
             $self->{'frame'},
             $self->{'orient'},
             $self->{'initn'},
             $self->{'init1'},
             $self->{'opt'},
             $self->{'zscore'},
             $self->{'bits'},
             $self->{'expect'},
            ) = (
                Bio::Parse::Format::FASTA::parse_frame($1),
                Bio::Parse::Format::FASTA::parse_orient($1),
                $2, $3, $4, $5, $6, $7,
            );
            next;
        }

        if ($line =~ /^
            (?:banded|trans.\s+)?
            (?:Smith-Waterman\s+score:\s*(\S+);)?    #sw score
            \s*(\S+)%\s+identity                     #percent identity
            \s+\((\S+)%\s+similar\)                  #percent similar
            \s+in\s+(\d+)                            #overlap length
            \s+(?:aa|nt)\s+overlap
            \s*
            (?:\s+\((\S+)\))                         #sequence ranges
            /xo) {

            $self->test_args(\$line,$2,$3);

            (
             $self->{'score'},
             $self->{'id_percent'},
             $self->{'similar_percent'},
             $self->{'overlap'},
             $self->{'ranges'},
            ) = (defined $1?$1:0,$2,$3,$4,defined $5?$5:'');
            next;
        }

        #should only get here for multiline descriptions
        if (defined $self->{'initn'}) {
            $self->warn("unknown field: $line");
            next;
        }

        #accumulate multiline descriptions (fasta... -L)
        $record .= ' ' . $line;
    }

    #now split out the description
    if ($record =~ /^
        >+
        (\S+)                      #id
        \s+
        (.*)                       #description
        \s+
        \(\s*(\d+)\s*(?:aa|nt)\)   #length
        \s*
        $/xo) {

        $self->test_args(\$record, $1, $2, $3);

        (
         $self->{'id'},
         $self->{'desc'},
         $self->{'length'},
        ) = (clean_identifier($1), strip_english_newlines($2), $3);

    } elsif ($record =~ /^>--/) {  #alternative alignment
        my $sib = $self->get_sibling(1);
        (
         $self->{'id'},
         $self->{'desc'},
         $self->{'length'},
         $self->{'alternative'},
        ) =
        (
         $sib->{'id'},
         $sib->{'desc'},
         $sib->{'length'},
         1,  #true
        );

    } else {
        $self->warn("unknown field: $record");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA3X::tfasta::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X::MATCH::ALN);

# tfast[axy]  pro x dna

sub query_orient { '+' }
sub sbjct_orient { $_[0]->get_summary->{'orient'} }

sub query_base { return 1 }
sub sbjct_base { return 3 }


###########################################################################
1;
