# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Handles: GCG FASTA   2.x
#
###########################################################################
package Bio::Parse::Format::FASTA2_GCG;

use Bio::Parse::Format::FASTA;
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

@ISA = qw(Bio::Parse::Format::FASTA);

@VERSIONS = (
             'GCG/2' => [
                 'FASTA',
             ],
            );

$NULL  = '^\s*$';

$ENTRY_START   = '^\(\S+\)\s+FASTA of:';
$ENTRY_END     = '! Output File:';

$HEADER_START  = $ENTRY_START;
$HEADER_END    = '^The best scores are:';

$RANK_START    = $HEADER_END;
$RANK_END      = $Bio::Parse::Format::FASTA::GCG_JUNK;

$TRAILER_START = '^! CPU time used:';
$TRAILER_END   = $ENTRY_END;

$MATCH_START   = '^\S+(?:\s+/rev)?\s*$';
$MATCH_END     = "(?:$MATCH_START|$TRAILER_START|$ENTRY_END)";

$SUM_START     = $MATCH_START;
$SUM_END       = '% identity in';

$ALN_START     = '^\s+\d+\s+';    #the ruler
$ALN_END       = '(?:^\S+(?:\s+/rev)?\s*$' . "|$MATCH_END)";


###########################################################################
#Parse one entry: generic for all FASTA family
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #Header lines
        if ($line =~ /$HEADER_START/o) {
            $scan->scan_until($HEADER_END);
            $self->push_record('HEADER',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Rank lines
        if ($line =~ /$RANK_START/o) {
            $scan->scan_until_inclusive($RANK_END);
            $self->push_record('RANK',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Hit lines
        if ($line =~ /$MATCH_START/o) {
            $scan->scan_skipping_until($MATCH_END, 1, 'MATCH');
            $self->push_record('MATCH',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Trailer lines
        if ($line =~ /$TRAILER_START/o) {
            $scan->scan_until_inclusive($TRAILER_END);
            $self->push_record('TRAILER',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #end of FASTA job
        next    if $line =~ /$ENTRY_END/o;

        #blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::HEADER);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'query'} = '';

    my ($ac, $id) = (undef, undef);

    while (defined ($line = $scan->read_line)) {

        #GCG fasta doesn't write explicit versioning

        #GCG filename and from/to giving query length
        if ($line =~ /^\s*
            \(\S+\)
            \s+
            FASTA of:\s+(\S+)    #file
            \s+
            from:\s+(\d+)        #from position
            \s+
            to:\s+(\d+)          #to position
            /xo) {

            $self->test_args(\$line, $1, $2);
            (
             $self->{'queryfile'},
             $self->{'length'},
            ) = ($1, ($3 > $2 ? $3-$2+1 : $2-$3+1));
            next;
        }

        #query identifier
        if ($line =~ /^\s*ID\s{3}(\S+);?/) {
            $self->test_args(\$line, $1);
            $id = $1;
            next;
        }

        #query accession
        if ($line =~ /^\s*AC\s{3}(\S+);?/) {
            $self->test_args(\$line, $1);
            $ac = $1;
            next;
        }

        #database size
        if ($line =~ /.*
            Sequences:\s+(\S+)       #sequence count (contains commas)
            \s+
            Symbols:\s+(\S+)         #symbol count (contains commas)
            /xo) {

            $self->test_args(\$line, $1,$2);

            (
             $self->{'sequences'},
             $self->{'residues'},
            ) = ($1, $2);

            $self->{'sequences'} =~ s/,//g;
            $self->{'residues'}  =~ s/,//g;

            next;
        }

        #ignore any other text

    }

    if (! defined $self->{'full_version'} ) {
        #can't determine version: hardwire one!
        $self->{'full_version'} = 'GCG FASTA 2';
        $self->{'version'}      = '2';
    }

    if (defined $id and defined $ac) {
        $self->{'query'} = "$ac:$id";
    } elsif (defined $id) {
        $self->{'query'} = "$id";
    } elsif (defined $ac) {
        $self->{'query'} = "$ac";
    } elsif (defined $self->{'queryfile'}) {
        $self->{'query'} = $self->{'queryfile'};
    } else {
        $self->{'query'} = 'query';
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'hit'} = [];

    #ranked search hits
    while (defined ($line = $scan->read_line)) {

        next    if $line =~ /$Bio::Parse::Format::FASTA2_GCG::RANK_START/o;
        next    if $line =~ /$Bio::Parse::Format::FASTA::GCG_JUNK/o;

        #first line
        if ($line =~ /^\s*
           ([^\s]+)                #id
           \s+
           Begin:\s+(\d+)          #start position
           \s+
           End:\s+(\d+)            #end position
           /xo) {

            #initialise
            my ($id,$desc,$init1,$initn,$opt,$z,$e) = ('','','','','','','');

            #warn "($1,$2,$3)\n";

            $self->test_args(\$line, $1);

            $id = clean_identifier($1);

            #read next line
            $line = $scan->read_line;

            #second line
            if ($line =~ /
                \s+
                (\d+)                   #init1
                \s+
                (\d+)                   #initn
                \s+
                (\d+)                   #opt
                \s+
                (\S+)                   #z-score
                \s+
                (\S+)                   #E(58765)
                \s*$
                $/xo) {

                #warn "($1,$2,$3,$4,$5)\n";

                $self->test_args(\$line, $1,$2,$3,$4,$5);

                ($init1,$initn,$opt,$z,$e) = ($1,$2,$3,$4,$5);

                $desc = $`;

                $desc =~ s/^\s*!\s*//;
            }

            push(@{$self->{'hit'}},
                 {
                  'id'     => $id,
                  'desc'   => $desc,
                  'initn'  => $initn,
                  'init1'  => $init1,
                  'opt'    => $opt,
                  'zscore' => $z,
                  'expect' => $e,
                 });

            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA2_GCG::NULL/o;

        #default
        $self->warn("unknown field: $line");
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::TRAILER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::TRAILER);


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #identifier lines
        if ($line =~ /$Bio::Parse::Format::FASTA2_GCG::SUM_START/o) {
            $scan->scan_until_inclusive($Bio::Parse::Format::FASTA2_GCG::SUM_END);
            $self->push_record('SUM',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #fragment hits: terminated by several possibilities
        if ($line =~ /$Bio::Parse::Format::FASTA2_GCG::ALN_START/o) {
            $scan->scan_until($Bio::Parse::Format::FASTA2_GCG::ALN_END);
            $self->push_record('ALN',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                    );
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA2_GCG::NULL/o;

        #ugly: skip queryfile or hit identifier lines
        next    if $line =~ /$Bio::Parse::Format::FASTA2_GCG::ALN_END/;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::MATCH::SUM;

use Bio::Parse::Strings qw(strip_english_newlines);
use Bio::Util::Regexp;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #skip first line (query filename)
    $scan->read_line;

    my ($id1, $id2, $ac, $de) = ('','','','');

    while (defined ($line = $scan->read_line)) {

        #hit database:identifier: appears as line 2
        if ($line =~ /^\s*(\S+:\S+);?\s*$/) {
            $id1 = $1;
            next;
        }

        #hit identifier
        if ($line =~ /^ID\s{3}(\S+);?.*;\s+(\d+)\s+\S+\./o) {
            $id2 = $1;
            $self->{'length'} = $2;
            next;
        }

        #hit accession
        if ($line =~ /^AC\s{3}(\S+);?/o) {
            $ac = $1;
            next;
        }

        #hit description
        if ($line =~ /^DE\s{3}(.*)/o) {
            $de .= strip_english_newlines($1);
            $de =~ s/\s*\.\s\.\s\.\s*$//;
            $de =~ s/\s*\.\.\.\*$//;
            next;
        }

        #skip other database entry lines
        next    if $line =~ /^(NI|DT|DE)\s{3}/o;

        #blank line or empty record: ignore
        next    if $line =~ /$Bio::Parse::Format::FASTA2_GCG::NULL/o;

        #scores
        if ($line =~ /^
            SCORES\s+
            init1\:\s*(\S+)        #init1
            \s*
            initn\:\s*(\S+)        #initn
            \s*
            opt\:\s*(\S+)          #opt
            \s*
            z-score\:\s*(\S+)      #z
            \s*
            E\(\)\:\s*(\S+)        #E
            \s*
            $/ixo) {

            $self->test_args(\$line,$1,$2,$3,$4,$5);

            (
             $self->{'init1'},
             $self->{'initn'},
             $self->{'opt'},
             $self->{'zscore'},
             $self->{'expect'},
            ) = ($1,$2,$3,$4,$5);

            $self->{'id'}   = ($id1 ne '' ? $id1 : $id2);
            $self->{'desc'} = $de;

            next;
        }

        if ($line =~ /^
            #smith-waterman in fasta2 and fasta3, maybe in gcg?
            (?:Smith-Waterman\s+score:\s*(\d+);)?    #sw score
            \s*($RX_Ureal)%                          #percent identity
            \s*identity\s+in\s+(\d+)                 #overlap length
            \s+(?:aa|nt|bp)\s+overlap
            \s*$/xo) {

            $self->test_args(\$line,$2,$3);

            (
             $self->{'score'},
             $self->{'id_percent'},
             $self->{'overlap'},
            ) = (defined $1?$1:0,$2,$3);

            next;
        }

        #default
        $self->warn("unknown field: $line");
    }

    $self;
}


###########################################################################
package Bio::Parse::Format::FASTA2_GCG::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA1::MATCH::ALN);


###########################################################################
1;
