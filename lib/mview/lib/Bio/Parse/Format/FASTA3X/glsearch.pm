# Copyright (C) 2013-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA3X::glsearch;

use Bio::Parse::Format::FASTA3X;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA3X::glsearch::HEADER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA3X::HEADER);

sub new {
    my $self = shift;
    $self = $self->SUPER::new(@_);
    #assume the query identifier is the same as the query filename
    if ($self->{query} eq '' and $self->{queryfile} ne '') {
        $self->{query} = $self->{queryfile};
    }
    return $self;
}


###########################################################################
package Bio::Parse::Format::FASTA3X::glsearch::RANK;

use Bio::Parse::Strings qw(clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #ranked search hits
    while (defined ($line = $scan->read_line)) {

        next    if $line =~ /$Bio::Parse::Format::FASTA3X::RANK_START/o;

        if($line =~ /^
           \s*
           ([^\s]+)                #id
           \s+
           (.*)                    #description: possibly empty
           \s+
           \(\s*(\S+)\)            #hit sequence length
           \s*
           (?:\[(\S)\])?           #frame
           \s*
           ($RX_Sint)              #n-w score
           \s+
           (\S+)                   #bits
           \s+
           (\S+)                   #E-value
           \s*
           $/xo) {

            $self->test_args(\$line, $1, $3, $5,$6,$7);

            push(@{$self->{'hit'}},
                 {
                  'id'      => clean_identifier($1),
                  'desc'    => $2,
                  #ignore $3
                  'frame'   => Bio::Parse::Format::FASTA::parse_frame($4),
                  'orient'  => Bio::Parse::Format::FASTA::parse_orient($4),
                  'score'   => $5,
                  'bits'    => $6,
                  'expect'  => $7,
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
package Bio::Parse::Format::FASTA3X::glsearch::TRAILER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA3X::glsearch::MATCH;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA3X::glsearch::MATCH::SUM;

use Bio::Parse::Strings qw(strip_english_newlines clean_identifier);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::MATCH::SUM);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);

    my $lines = $scan->read_until_inclusive('^\s*global/local');

    if ($lines =~ /^
        >*
        (\S+)                             #id
        \s+
        (.*)                              #description: possibly empty
        \s+
        \(\s*(\d+)\s*(?:aa|nt)\)          #length
        \s*
        (rev-comp)?                       #frame
        \s+
        n-w\s+opt:\s*($RX_Sint)           #n-w opt
        \s+
        Z-score:\s*(\S+)                  #Z-score
        \s+
        bits:\s*(\S+)                     #bits
        \s+
        E\((?:\d+)?\):\s*(\S+)            #expect
        \s+
        global\/local\s+score:\s*($RX_Sint);   #n-w score (same as opt?)
        \s+
        (\S+)%\s+identity                 #percent identity
        \s+
        \((\S+)%\s+(?:similar|ungapped)\) #percent similarity
        \s+in\s+
        (\d+)\s*(?:aa|nt)\s+overlap       #alignment overlap
        \s+
        \((\S+)\)                         #alignment ranges
        \s*
        $/xso) {

        $self->test_args(\$lines, $1, $3, $5,$6,$7,$8,$9,$10,$11,$12,$13);

        (
         $self->{'id'},
         $self->{'desc'},
         $self->{'length'},
         $self->{'frame'},
         $self->{'orient'},
         $self->{'opt'},
         $self->{'zscore'},
         $self->{'bits'},
         $self->{'expect'},
         $self->{'score'},
         $self->{'id_percent'},
         $self->{'sim_percent'},
         $self->{'overlap'},
         $self->{'ranges'},
        ) = (
            clean_identifier($1),
            strip_english_newlines($2),
            $3,
            Bio::Parse::Format::FASTA::parse_frame($4),
            Bio::Parse::Format::FASTA::parse_orient($4),
            $5, $6, $7, $8, $9, $10, $11, $12, $13,
        );

    } elsif ($lines =~ /^>--/) {  #alternative alignment
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
        $self->warn("unknown field: $lines");
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA3X::glsearch::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA3X::MATCH::ALN);

sub query_orient { $_[0]->get_summary->{'orient'} }
sub sbjct_orient { '+' }


###########################################################################
1;
