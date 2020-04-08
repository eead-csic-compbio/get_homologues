# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Handles: fasta, tfastx  2.x
#
###########################################################################
package Bio::Parse::Format::FASTA2;

use Bio::Parse::Format::FASTA;
use Bio::Parse::Format::FASTA2::align;
use Bio::Parse::Format::FASTA2::lalign;

use strict;

use vars qw(
            @ISA

            @VERSIONS

            $NULL

            $ENTRY_START
            $ENTRY_END

            $HEADER_START
            $HEADER_END

            $RANK_START
            $RANK_END

            $TRAILER_START
            $TRAILER_END

            $MATCH_START
            $MATCH_END

            $SUM_START
            $SUM_END

            $ALN_START
            $ALN_END
);

@ISA   = qw(Bio::Parse::Format::FASTA);

@VERSIONS = (
             '2' => [
                     'FASTA',
                     'TFASTX',
                     'ALIGN',
                     'LALIGN',
                     'SSEARCH',
                    ],
            );

$NULL  = '^\s*$';

$ENTRY_START   = '(?:'
    . '^\s*FASTA searches a protein or DNA sequence data bank'
    . '|'
    . '^\s*TFASTX translates and searches a DNA sequence data bank'
    . '|'
    . '^\s*\S+\s*[,:]\s+\d+\s+(?:aa|nt)'
    . '|'
    . 'SSEARCH searches a sequence database'
    . '|'
    . $Bio::Parse::Format::FASTA2::align::ALIGN_START
    . '|'
    . $Bio::Parse::Format::FASTA2::lalign::ALIGN_START
    . ')';
$ENTRY_END     = '(?:'
    . 'Library scan:'
    . '|'
    . $Bio::Parse::Format::FASTA2::align::ALIGN_END
    . '|'
    . $Bio::Parse::Format::FASTA2::lalign::ALIGN_END
    . ')';

$HEADER_START  = $ENTRY_START;
$HEADER_END    = '^The best scores are:';

$RANK_START    = $HEADER_END;
$RANK_END      = $NULL;

$TRAILER_START = $ENTRY_END;
$TRAILER_END   = $ENTRY_END;

$MATCH_START   = '^>*\S{7}.*\(\d+ (?:aa|nt)\)';
$MATCH_END     = "(?:$MATCH_START|$ENTRY_END)";

$SUM_START     = $MATCH_START;
$SUM_END       = $NULL;

$ALN_START     = '^(?:\s+\d+\s+|\s+$)';  #the ruler
$ALN_END       = $MATCH_END;

#Parse one entry
sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA2::HEADER;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Format::FASTA::HEADER);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'query'}     = '';
    $self->{'queryfile'} = '';

    while (defined ($line = $scan->read_line)) {

        if ($line =~ /^\s*(version\s+(\S+).*)/) {
            $self->{'full_version'} = $1;
            $self->{'version'}      = $2;
            next;
        }

        if ($line =~ /^\s*v((\d+(?:\.\d+)?)\s*\S+.*)/) {
            $self->{'full_version'} = $1;
            $self->{'version'}      = $2;
            next;
        }

        if ($line =~ /^\s*>?(\S+).*\:\s+(\d+)\s+(?:aa|nt)/) {
            if ($self->{'queryfile'} eq '') {
                $self->{'queryfile'} = $1;
                $self->{'queryfile'} =~ s/,$//;
            } else {
                $self->test_args(\$line, $1);
                $self->{'query'} = clean_identifier($1);
            }
            $self->{'length'} = $2;
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
        $self->{'full_version'} = 'looks like FASTA 2';
        $self->{'version'}      = '2.1';
        $self->warn("unknown FASTA version - using $self->{'version'}");    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA2::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA1::MATCH::ALN);


###########################################################################
1;
