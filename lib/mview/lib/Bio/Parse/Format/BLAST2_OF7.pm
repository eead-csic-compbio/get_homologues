# Copyright (C) 2015-2020 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Base classes for NCBI BLAST2 family.
#
# Handles: BLAST+ -outfmt 7
#
# BLAST (NCBI version 2) iterated searching uses 3 main record types:
#
#   HEADER        the header text
#   SEARCH {*}    passes of the search engine
#     RANK        the list of ordered high scoring hits
#     MATCH {*}   the set of alignments for a given hit
#       SUM         the summary lines for each hit
#       ALN {*}     each aligned fragment: score + alignment
#
###########################################################################
package Bio::Parse::Format::BLAST2_OF7;

use Bio::Parse::Format::BLAST;
use Bio::Parse::Format::BLAST2;

use strict;

use vars qw(@ISA

            @VERSIONS

            $ENTRY_START
            $ENTRY_END

            $HEADER_START
            $HEADER_END
           );

@ISA = qw(Bio::Parse::Format::BLAST);

@VERSIONS = (
             '2of7' => [
                     'BLASTP',
                     'BLASTN',
                     'BLASTX',
                     'TBLASTN',
                     'TBLASTX',
                     'PSIBLAST',
                     'PHIBLASTP',
                   ],
            );

# BLAST -outfmt 7 HEADER format is:
#
# Upto blast version 2.4.0+
#
#  # PSIBLAST 2.2.28+
#  # Iteration: 1                  [note: psiblast only]
#  # Query: test
#  # Database: mito.1000.aa
#  # Fields: query id, subject id, % identity, alignment length, ...
#  # 65 hits found
# <data>\t<data>\t...
# <data>\t<data>\t...
# ...

# From blast version 2.5.0+
#
#  # PSIBLAST 2.5.0+
#  # Iteration: 1                  [note: psiblast only]
#  # Query: test
#  # Database: mito.1000.aa
#  # Fields: query acc., subject acc., % identity, alignment length, ...
#  # 65 hits found
# <data>\t<data>\t...
# <data>\t<data>\t...
# ...

# From blast version 2.6.0+
#
#  # PSIBLAST 2.6.0+
#  # Iteration: 1                  [note: psiblast only]
#  # Query: test
#  # Database: mito.1000.aa
#  # Fields: query acc.ver, subject acc.ver, % identity, alignment length, ...
#  # 65 hits found
# <data>\t<data>\t...
# <data>\t<data>\t...
# ...

# followed by
#
# default case: single search terminates with:
#  # BLAST processed 1 queries
#  <eof>

# psiblast searches repeat then terminate with:
#
#  <blank>
#  Search has CONVERGED!
#  # BLAST processed 2 queries
#  <eof>

my $NULL = '^\s*$';

my $PROGRAMS = "(?:" . join("|", @{$VERSIONS[1]}) . ")";

$ENTRY_START     = "^\# $PROGRAMS";
$ENTRY_END       = '^\# .*processed';

$HEADER_START    = $ENTRY_START;
$HEADER_END      = '^[^\#]';

my $SEARCH_START = "^[^\#][^\t]*\t";
my $SEARCH_END   = "^(?:$NULL|$HEADER_START|$ENTRY_END)";

my $FIELD_SKIP   = '-';
my $FIELD_MAP    = {  #blastp -help for -outfmt 7 fields
    #HEADER comment                   blast CL specifier   #CL 'std' specifier
    #-----------------------------------------------------------------------
    'query id'                        => 'qseqid',         #upto 2.4.0+
    'query gi'                        => 'qgi',
    'query acc.'                      => 'qacc',           #from 2.5.0+
    'query acc.ver'                   => 'qaccver',        #from 2.6.0+
    'query length'                    => 'qlen',
    'subject id'                      => 'sseqid',         #upto 2.4.0+
    'subject ids'                     => 'sallseqid',
    'subject gi'                      => 'sgi',
    'subject gis'                     => 'sallgi',
    'subject acc.'                    => 'sacc',           #from 2.5.0+
    'subject acc.ver'                 => 'saccver',        #from 2.6.0+
    'subject accs.'                   => 'sallacc',
    'subject length'                  => 'slen',
    'q. start'                        => 'qstart',         #yes
    'q. end'                          => 'qend',           #yes
    's. start'                        => 'sstart',         #yes
    's. end'                          => 'send',           #yes
    'query seq'                       => 'qseq',
    'subject seq'                     => 'sseq',
    'evalue'                          => 'evalue',         #yes
    'bit score'                       => 'bitscore',       #yes
    'score'                           => 'score',
    'alignment length'                => 'length',         #yes
    '% identity'                      => 'pident',         #yes
    'identical'                       => 'nident',
    'mismatches'                      => 'mismatch',       #yes
    'positives'                       => 'positive',
    'gap opens'                       => 'gapopen',        #yes
    'gaps'                            => 'gaps',
    '% positives'                     => 'ppos',
    'query/sbjct frames'              => 'frames',
    'query frame'                     => 'qframe',
    'sbjct frame'                     => 'sframe',
    'BTOP'                            => 'btop',
    'subject tax id'                  => 'staxid',         #from 2.4.0+
    'subject sci name'                => 'ssciname',       #from 2.4.0+
    'subject com names'               => 'scomname',       #from 2.4.0+
    'subject blast name'              => 'sblastname',     #from 2.4.0+
    'subject super kingdom'           => 'sskingdom',      #from 2.4.0+
    'subject tax ids'                 => 'staxids',
    'subject sci names'               => 'sscinames',
    'subject com names'               => 'scomnames',  #blast: header reused
    'subject blast names'             => 'sblastnames',
    'subject super kingdoms'          => 'sskingdoms',
    'subject title'                   => 'stitle',
    'subject titles'                  => 'salltitles',
    'subject strand'                  => 'sstrand',
    '% subject coverage'              => 'qcovs',
    '% hsp coverage'                  => 'qcovhsp',
    'query coverage per uniq subject' => 'qcovus',         #from 2.3.0+
};

#required fields: blast option name to mview internal key
my $MAP_RANK = {
    'sseqid'      => 'id',    #upto 2.4.0+
    'sacc'        => 'id',    #from 2.5.0+
    'saccver'     => 'id',    #from 2.6.0+
    'evalue'      => 'expect',
    'bitscore'    => 'bits',
    'stitle'      => 'summary',
    'salltitles'  => 'summary',
};

my $MAP_SUM = {
    'sseqid'      => 'id',    #upto 2.4.0+
    'sacc'        => 'id',    #from 2.5.0+
    'saccver'     => 'id',    #from 2.6.0+
    'length'      => 'length',
    'stitle'      => 'desc',
    'salltitles'  => 'desc',
};

my $MAP_ALN = {
    'evalue'      => 'expect',
    'bitscore'    => 'bits',
    'qseq'        => 'query',
    'qstart'      => 'query_start',
    'qend'        => 'query_stop',
    'sseq'        => 'sbjct',
    'sstart'      => 'sbjct_start',
    'send'        => 'sbjct_stop',
    'pident'      => 'id_percent',
};

#pass-through fields: blast option name
my $MAP_EXTRA = {
    'staxid'      => 1,
    'ssciname'    => 1,
    'scomname'    => 1,
    'sblastname'  => 1,
    'sskingdom'   => 1,
    'staxids'     => 1,
    'sscinames'   => 1,
    'scomnames'   => 1,
    'sblastnames' => 1,
    'sskingdoms'  => 1,
};

#Given a string in 'line' of tab-separated fields named as in 'blastopts',
#extract those in 'wanted' storing each such key/value into 'hash'; returns
#number of fields read or -1 on error.
sub get_fields {
    my ($hash, $wanted, $line, $fields, $counts, $debug) = (@_, 0);

    my @list = split("\t", $line);

    if ($debug) {
        #warn "GF: [@$fields] => [@{[keys %$wanted]}] => [@{[values %$wanted]}]\n";
        #warn "GF: [$line]\n";
    }
    return -1  if scalar @list != scalar @$fields;

    my $c = 0;
    for (my $i=0; $i < @$fields; $i++) {
        my $blastopt = $fields->[$i];
        my $val = $list[$i];

        $val =~ s/^\s+|\s+$//g;  #blast 2.2.28+, strip extra space

        #ignore empty or missing data
        next  if $val eq 'N/A' or $val eq 'n/a' or $val eq '';

        $counts->[$i]++  if defined $counts;

        if (exists $wanted->{$blastopt}) {  #required key/value
            my $mview_attr = $wanted->{$blastopt};
            $hash->{$mview_attr} = $val;
            warn "GF: [$blastopt] => [$mview_attr] => [$val]\n"  if $debug;
            $c++;  #count required fields
        } elsif (exists $MAP_EXTRA->{$blastopt}) {  #optional key/value
            $hash->{'extra'}->{$blastopt} = $val;
            warn "GF: [extra] => [$blastopt] => [$val]\n"  if $debug;
        }
    }
    #warn "GF: read $c fields\n"  if $debug;
    return $c;
}

#strip leading identifier from summary string
sub strip_id {
    my ($id, $s) = @_;
    $s = substr($s, length($id))  if index($s, $id) == 0;
    $s =~ s/^\s+|\s+$//g;  #strip extra space
    $s;
}

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #blank line or empty record: ignore
        next  if $line =~ /$NULL/o;

        #HEADER block
        if ($line =~ /$HEADER_START/o) {
            $scan->scan_until($HEADER_END);
                $self->push_record('HEADER',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #SEARCH block
        if ($line =~ /$SEARCH_START/o) {
            $scan->scan_until($SEARCH_END);
                $self->push_record('SEARCH',
                                   $scan->get_block_start(),
                                   $scan->get_block_stop(),
                    );
            next;
        }

        #stop at psiblast convergence message
        last  if $line =~ /^Search has CONVERGED/;

        #stop before terminal comment
        last  if $line =~ /$ENTRY_END/;

        #default
        $self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #BLAST
    $self->{'full_version'} = '';
    $self->{'version'}      = '';
    $self->{'query'}        = '';
    $self->{'summary'}      = '';

    #BLAST2_OF7
    $self->{'fields'}       = [];
    $self->{'counts'}       = [];

    while (defined ($line = $scan->read_line(1))) {

        #blast version info
        if ($line =~ /^# ($PROGRAMS\s+(\S+))/o) {
            $self->test_args(\$line, $1, $2);
            (
             $self->{'full_version'},
             $self->{'version'},
            ) = ($1, $2);
            next;
        }

        if ($line =~ /^# Query:\s+(.*)/o) {
            $self->{'query'} = $1;
            $self->{'summary'} = '';  #never set
            next;
        }

        if ($line =~ /^# Fields:\s+(.*)/o) {
            $self->init_fields($1);
            $self->extract_fields;
            next;
        }

        next  if $line =~ /^# Iteration:\s+(\d+)/o;
        next  if $line =~ /^# Database:/o;
        next  if $line =~ /^# \d+ hits found/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}

sub init_fields {
    my $self = shift;
    foreach my $f (split(/,\s+/, $_[0])) {
        my $blastopt = $FIELD_SKIP;
        $blastopt = $FIELD_MAP->{$f}  if exists $FIELD_MAP->{$f};
        push @{$self->{'fields'}}, $blastopt;
        push @{$self->{'counts'}}, 0;
    }
}

sub extract_fields {
    my $self = shift;
    foreach my $blastopt (@{$self->{'fields'}}) {
        $self->{'extra'}->{$blastopt} = ''  if exists $MAP_EXTRA->{$blastopt};
    }
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> %s\n", 'version',      $self->{'version'};
    $s .= sprintf "$x%20s -> %s\n", 'full_version', $self->{'full_version'};
    $s .= sprintf "$x%20s -> %s\n", 'query',        $self->{'query'};
    $s .= sprintf "$x%20s -> %s\n", 'summary',      $self->{'summary'};
    $s .= sprintf "$x%20s -> %s\n", 'fields',       $self->fmt($self->{'fields'});
    $s .= sprintf "$x%20s -> %s (initial)\n", 'counts', $self->fmt($self->{'counts'});
    return $s;
}


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::SEARCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH);

sub new {
    my $self = new Bio::Parse::Record(@_);

    my ($type, $parent, $text, $start, $stop) = @_;

    #create SEARCH::RANK
    $self->push_record('RANK', $start, $stop);

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::SEARCH::RANK;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::RANK);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    #column headers
    $self->{'header'} = '';

    #ranked search hits
    $self->{'hit'}    = [];

    #fetch fields
    my $fields = $self->get_parent(2)->get_record('HEADER')->{'fields'};
    my $counts = $self->get_parent(2)->get_record('HEADER')->{'counts'};

    $self->init_field_counts($counts);

    #warn "[@{[scalar @$fields]}] [@{$fields}]\n";
    #warn "[@{[scalar @$counts]}] [@{$counts}]\n";

    #accumulate text block data for same hit id on successive lines
    my ($mid, $mstart, $mstop) = ('', 0, 0);

    my $parent = $self->get_parent(1);

    while (defined ($line = $scan->read_line(1))) {
        #warn "[$line]\n";

        my $tmp = {};

        if ($line =~ /$SEARCH_START/o) {

            #BLAST2
            $tmp->{'id'}      = '';
            $tmp->{'bits'}    = '';
            $tmp->{'expect'}  = '';
            $tmp->{'n'}       = '0';  #default
            $tmp->{'summary'} = '';

            #extract relevant fields
            my $c = Bio::Parse::Format::BLAST2_OF7::get_fields($tmp,
                                                               $MAP_RANK,
                                                               $line,
                                                               $fields,
                                                               $counts);

            if ($c < 0) {
                $self->die("field count mismatch (expect @{[scalar @$fields]}, got $c)\n");
            }

            if ($tmp->{'id'} eq '') {
                $self->die("blast column 'sseqid/sacc' is needed to identify hits");
            }

            $tmp->{'summary'} =
                Bio::Parse::Format::BLAST2_OF7::strip_id($tmp->{'id'},
                                                         $tmp->{'summary'});
            $tmp->{'id'} = clean_identifier($tmp->{'id'});

            #id same as last line: extend SEARCH::MATCH block
            if ($mid eq $tmp->{'id'}) {
                $parent->pop_record();
                $mstop = $scan->get_line_stop();
                $parent->push_record('MATCH', $mstart, $mstop);
                #expect the first hit supplied by blast to be the highest
                #scoring, but test just in case:
                if ($tmp->{'bits'} > $self->{'hit'}->[-1]->{'bits'}) {
                    #warn "later hit has greater bits score\n";
                    pop  @{$self->{'hit'}};
                    push @{$self->{'hit'}}, $tmp;
                }
                next;
            }

            #new id: create SEARCH::MATCH
            ($mid, $mstart, $mstop) = ($tmp->{'id'},
                                       $scan->get_line_start(),
                                       $scan->get_line_stop());
            $parent->push_record('MATCH', $mstart, $mstop);

            push @{$self->{'hit'}}, $tmp;

            next;
        }

        #blank line or empty record: ignore
        next  if $line =~ /$NULL/o;

        #psiblast convergence message: ignore
        next  if $line =~ /^Search has CONVERGED/o;

        #default
        $self->warn("unknown field: $line");
    }
    #warn "RANK(counts): [@{[scalar @$counts]}] [@{$counts}]\n";

    $self;#->examine;
}

sub init_field_counts {
    my ($self, $counts) = @_;
    for (my $i=0; $i < @$counts; $i++) {
        $counts->[$i] = 0;
    }
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> '%s'\n", 'header', $self->{'header'};
    foreach my $hit (@{$self->{'hit'}}) {
        foreach my $field (sort keys %$hit) {
            my $val = $hit->{$field};
            $val = $self->fmt_hash($val)  if $field eq 'extra';
            $s .= sprintf "$x%20s -> %s\n", $field, $val;
        }
    }
    return $s;
}


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST::MATCH);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    my ($type, $parent, $text, $start, $stop) = @_;

    #create SEARCH::MATCH::SUM
    $self->push_record('SUM', $start, $stop);

    while (defined ($line = $scan->read_line(1))) {
        #warn "[$line]\n";

        if ($line =~ /$SEARCH_START/o) {
            #create SEARCH::MATCH::ALN
            $self->push_record('ALN',
                               $scan->get_line_start(),
                               $scan->get_line_stop());
            next;
        }

        #default
        $self->warn("unknown field: $line");
    }

    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::SUM;

use Bio::Parse::Strings qw(clean_identifier);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::SUM);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);

    #BLAST
    $self->{'id'}     = '';
    $self->{'desc'}   = '';
    $self->{'length'} = '';

    $self->extract_fields($MAP_SUM, $scan->read_line(1));

    $self->{'desc'} =
        Bio::Parse::Format::BLAST2_OF7::strip_id($self->{'id'},
                                                 $self->{'desc'});

    $self->{'id'} = clean_identifier($self->{'id'});

    $self;#->examine;
}

sub extract_fields {
    my ($self, $want, $line) = @_;
    my $record = $self->get_parent(3)->get_record('HEADER');
    Bio::Parse::Format::BLAST2_OF7::get_fields($self, $want, $line,
                                               $record->{'fields'}, undef);
}


###########################################################################
package Bio::Parse::Format::BLAST2_OF7::SEARCH::MATCH::ALN;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::BLAST2::SEARCH::MATCH::ALN);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);

    #BLAST
    $self->{'query'}        = '';
    $self->{'align'}        = '';
    $self->{'sbjct'}        = '';
    $self->{'query_orient'} = '';
    $self->{'query_start'}  = '';
    $self->{'query_stop'}   = '';
    $self->{'sbjct_orient'} = '';
    $self->{'sbjct_start'}  = '';
    $self->{'sbjct_stop'}   = '';

    #BLAST2
    $self->{'bits'}         = '';
    $self->{'score'}        = '';
    $self->{'n'}            = '1';  #default
    $self->{'expect'}       = '';
    $self->{'id_fraction'}  = '';
    $self->{'id_percent'}   = '';
    $self->{'pos_fraction'} = '';
    $self->{'pos_percent'}  = '';
    $self->{'gap_fraction'} = '';
    $self->{'gap_percent'}  = '';

    $self->extract_fields($MAP_ALN, $scan->read_line(1));

    #use sequence numbering to get orientations
    $self->{'query_orient'} =
        $self->{'query_start'} > $self->{'query_stop'} ? '-' : '+';
    $self->{'sbjct_orient'} =
        $self->{'sbjct_start'} > $self->{'sbjct_stop'} ? '-' : '+';

    $self;#->examine;
}

sub extract_fields {
    my ($self, $want, $line) = @_;
    my $record = $self->get_parent(3)->get_record('HEADER');
    Bio::Parse::Format::BLAST2_OF7::get_fields($self, $want, $line,
                                               $record->{'fields'}, undef);
}


###########################################################################
1;
