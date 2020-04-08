# Copyright (C) 2013-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::FASTA3X::fastm;

use Bio::Parse::Format::FASTA3X;
use Bio::Parse::Format::FASTA3X::fasta;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X);

sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA3X::fastm::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X::fasta::HEADER);


###########################################################################
package Bio::Parse::Format::FASTA3X::fastm::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #ranked search hits
    while (defined ($line = $scan->read_line)) {

        next    if $line =~ /$Bio::Parse::Format::FASTA3X::RANK_START/o;

        #fasta3X behaviour
        if ($line =~ /^
            \s*
            (\S+)                #id
            \s+
            (.*)                 #description (may be empty)
            \s+
            (?:\[[^]]+\])?       #don't know - reported by rls@ebi.ac.uk
            \s*
            \(\s*(\d+)\)         #aa
            \s*
            (?:\[(\S)\])?        #frame
            \s+
            (\d+)                #initn
            \s+
            (\d+)                #init1
            \s+
            (\S+)                #bits
            \s+
            (\S+)                #E(205044)
            \s+
            (\S+)                #sn
            \s+
            (\S+)?               #sl
            \s*
            $/xo) {

            $self->test_args(\$line, $1, $3, $5,$6,$7); #not $2,$4

            push(@{$self->{'hit'}},
                 {
                  'id'     => clean_identifier($1),
                  'desc'   => $2,
                  'length' => $3,
                  'frame'  => Bio::Parse::Format::FASTA::parse_frame($4),
                  'initn'  => $5,
                  'init1'  => $6,
                  'bits'   => $7,
                  'expect' => $8,
                  'sn'     => $9,
                  'sl'     => $10,
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
package Bio::Parse::Format::FASTA3X::fastm::TRAILER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X::fasta::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA3X::fastm::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X::fasta::MATCH);


###########################################################################
package Bio::Parse::Format::FASTA3X::fastm::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA3X::fasta::MATCH::SUM);


###########################################################################
package Bio::Parse::Format::FASTA3X::fastm::MATCH::ALN;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA3X::MATCH::ALN);

#Although this a S-W alignment, the sequence range numbers look more
#like those of a global/global alignment: query sequence range is complete;
#sbjct depends on it
sub get_start_stop {
    my ($self, $tgt, $orient, $info, $base) = @_;

    my ($start_info, $stop_info, $fs1, $fs2) = @$info;

    my ($r1) = @{$start_info};
    my ($r2) = @{$stop_info};

    warn "$tgt(o): $orient,$base start/stop: $r1,$r2\n\n" if $Bio::Parse::Format::FASTA::DEBUG;
    return ($r1, $r2);
}


###########################################################################
1;
