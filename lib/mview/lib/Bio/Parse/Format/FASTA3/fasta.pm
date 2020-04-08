# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA3::fasta;

use Bio::Parse::Format::FASTA3;
use strict;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA3::fasta::HEADER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA3::HEADER);


###########################################################################
package Bio::Parse::Format::FASTA3::fasta::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #ranked search hits
    while (defined ($line = $scan->read_line)) {

        next    if $line =~ /$Bio::Parse::Format::FASTA3::RANK_START/o;

        #pre-fasta3.3 behaviour
        if ($line =~ /^
            \s*
            (\S+)                #id
            \s+
            (.*)                 #description (may be empty)
            \s+
            (?:\[[^]]+\])?       #don't know - reported by rls@ebi.ac.uk
            \s*
            \(\s*(\d+)\)         #aa
            \s+
            (?:\[(\S)\])?        #orientation
            \s*
            (\d+)                #initn
            \s+
            (\d+)                #init1
            \s+
            (\d+)                #opt
            \s+
            (\S+)                #z-score
            \s+
            (\S+)                #E(205044)
            \s*
            $/xo) {

            $self->test_args(\$line, $1, $3, $5,$6, $7,$8,$9); #not $2,$4

            push(@{$self->{'hit'}},
                 {
                  'id'     => clean_identifier($1),
                  'desc'   => $2,
                  'length' => $3,
                  'frame'  => Bio::Parse::Format::FASTA::parse_frame($4),
                  'orient' => Bio::Parse::Format::FASTA::parse_orient($4),
                  'initn'  => $5,
                  'init1'  => $6,
                  'opt'    => $7,
                  'zscore' => $8,
                  'bits'   => 0,
                  'expect' => $9,
                 });
            next;
        }

        #fasta3.3 behaviour
        if ($line =~ /^
            \s*
            (\S+)                #id
            \s+
            (.*)                 #description (may be empty)
            \s+
            (?:\[[^]]+\])?       #don't know - reported by rls@ebi.ac.uk
            \s*
            \(\s*(\d+)\)         #aa
            \s+
            (?:\[(\S)\])?        #orientation
            \s*
            (\d+)                #opt
            \s+
            (\S+)                #bits
            \s+
            (\S+)                #E(205044)
            \s*
            $/xo) {

            $self->test_args(\$line, $1, $3, $5,$6,$7); #not $2,$4

            push(@{$self->{'hit'}},
                 {
                  'id'     => clean_identifier($1),
                  'desc'   => $2,
                  'length' => $3,
                  'frame'  => Bio::Parse::Format::FASTA::parse_frame($4),
                  'orient' => Bio::Parse::Format::FASTA::parse_orient($4),
                  'initn'  => 0,
                  'init1'  => 0,
                  'opt'    => $5,
                  'zscore' => 0,
                  'bits'   => $6,
                  'expect' => $7,
                 });
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA3::NULL/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA3::fasta::TRAILER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA3::fasta::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA3::fasta::MATCH::SUM;

use Bio::Parse::Strings qw(strip_english_newlines clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::MATCH::SUM);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    my $record = '';

    while (defined ($line = $scan->read_line(1))) {

        if ($line =~ /^>+/) {
            $record = $line;       #process this later
            next;
        }

        #pre-fasta3.3 behaviour
        if ($line =~ /^
            (rev-comp)?            #frame
            \s*
            initn\:\s*(\S+)        #initn
            \s*
            init1\:\s*(\S+)        #init1
            \s*
            opt\:\s*(\S+)          #opt
            \s*
            Z-score\:\s*(\S+)      #z
            \s*
            expect\(\)\s*(\S+)     #E
            \s*
            $/xo) {

            $self->test_args(\$line,$2,$3,$4,$5,$6); #not $1

            (
             $self->{'frame'},
             $self->{'orient'},
             $self->{'initn'},
             $self->{'init1'},
             $self->{'opt'},
             $self->{'zscore'},
             $self->{'expect'},
            ) = (
                Bio::Parse::Format::FASTA::parse_frame($1),
                Bio::Parse::Format::FASTA::parse_orient($1),
                $2, $3, $4, $5, $6,
            );
            next;
        }

        #fasta3.3 behaviour
        if ($line =~ /^
            (rev-comp)?            #frame
            \s*
            initn\:\s*(\S+)        #initn
            \s*
            init1\:\s*(\S+)        #init1
            \s*
            opt\:\s*(\S+)          #opt
            \s*
            Z-score\:\s*(\S+)      #z
            \s*
            bits\:\s*(\S+)         #bits
            \s*
            E\((?:\d+)?\):\s*(\S+) #E; fasta36 adds parenthesized integer
            \s*
            /xo) {

            $self->test_args(\$line,$2,$3,$4,$5,$6,$7); #not $1

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
            (?:(?:banded\s+)?Smith-Waterman\s+score:\s*(\d+);)?  #sw score
            \s*($RX_Ureal)%\s*identity               #percent identity
            (?:\s+\(($RX_Ureal)%\s+(?:ungapped|similar)\))?  #percent similar
            \s+in\s+(\d+)                            #overlap length
            \s+(?:aa|nt)\s+overlap
            (?:\s+\((\S+)\))?                        #sequence ranges
            /xo) {

            $self->test_args(\$line,$2,$4);

            (
             $self->{'score'},
             $self->{'id_percent'},
             $self->{'id_percent_similar'},
             $self->{'overlap'},
             $self->{'ranges'},
            ) = (defined $1?$1:0,$2,defined $3?$3:'',$4,defined $5?$5:'');
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA3::NULL/o;

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
        >*
        (\S+)                       #id
        \s+
        (.*)                        #description (may be empty)
        \s+
        \(\s*(\d+)\s*(?:aa|nt)\)    #length
        /xo) {

        $self->test_args(\$record, $1, $3); #not $2

        (
         $self->{'id'},
         $self->{'desc'},
         $self->{'length'},
         ) = (clean_identifier($1), strip_english_newlines($2), $3);

    } else {
        $self->warn("unknown field: $record");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA3::fasta::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA3::MATCH::ALN);

sub query_orient { $_[0]->get_summary->{'orient'} }
sub sbjct_orient { '+' }


###########################################################################
1;
