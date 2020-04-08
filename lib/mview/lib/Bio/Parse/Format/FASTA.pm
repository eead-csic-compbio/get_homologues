# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Base classes for FASTA family.
#
# FASTA parsing consists of 4 main record types:
#
#   HEADER       the header text (not much interesting)
#   RANK         the list of ordered hits and initn/init1/opt (scores)
#   MATCH        each hit's alignment with the query
#   TRAILER      trailing run-time data
#
# MATCH is further subdivided into:
#   SUM          the summary lines for each hit: name, description, scores
#   ALN          the aligned fragments
#
###########################################################################
package Bio::Parse::Format::FASTA;

use vars qw(@ISA $GCG_JUNK);
use strict;

@ISA = qw(Bio::Parse::Record);

BEGIN { $GCG_JUNK = '(?:^\.\.|^\\\\)' }

use Bio::Parse::Format::FASTA1;
use Bio::Parse::Format::FASTA2;
use Bio::Parse::Format::FASTA3;
use Bio::Parse::Format::FASTA3X;
use Bio::Parse::Format::FASTA2_GCG;

my %VERSIONS = (
                @Bio::Parse::Format::FASTA1::VERSIONS,
                @Bio::Parse::Format::FASTA2::VERSIONS,
                @Bio::Parse::Format::FASTA3::VERSIONS,
                @Bio::Parse::Format::FASTA3X::VERSIONS,
                @Bio::Parse::Format::FASTA2_GCG::VERSIONS,
               );

my $NULL        = '^\s*$';

my $FASTA_START = '^\s*\S+\s*[,:]\s+\d+\s+(?:aa|nt)';

my $ENTRY_START = "(?:"
    . $Bio::Parse::Format::FASTA1::ENTRY_START
    . "|"
    . $Bio::Parse::Format::FASTA2::ENTRY_START
    . "|"
    . $Bio::Parse::Format::FASTA3::ENTRY_START
    . "|"
    . $Bio::Parse::Format::FASTA3X::ENTRY_START
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::ENTRY_START
    . ")";
my $ENTRY_END = "(?:"
    . $Bio::Parse::Format::FASTA1::ENTRY_END
    . "|"
    . $Bio::Parse::Format::FASTA2::ENTRY_END
    . "|"
    . $Bio::Parse::Format::FASTA3::ENTRY_END
    . "|"
    . $Bio::Parse::Format::FASTA3X::ENTRY_END
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::ENTRY_END
    . ")";

my $HEADER_START = "(?:"
    . $Bio::Parse::Format::FASTA1::HEADER_START
    . "|"
    . $Bio::Parse::Format::FASTA2::HEADER_START
    . "|"
    . $Bio::Parse::Format::FASTA3::HEADER_START
    . "|"
    . $Bio::Parse::Format::FASTA3X::HEADER_START
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::HEADER_START
    . ")";
my $HEADER_END = "(?:"
    . $Bio::Parse::Format::FASTA1::HEADER_END
    . "|"
    . $Bio::Parse::Format::FASTA2::HEADER_END
    . "|"
    . $Bio::Parse::Format::FASTA3::HEADER_END
    . "|"
    . $Bio::Parse::Format::FASTA3X::HEADER_END
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::HEADER_END
    . ")";

my $RANK_START = "(?:"
    . $Bio::Parse::Format::FASTA1::RANK_START
    . "|"
    . $Bio::Parse::Format::FASTA2::RANK_START
    . "|"
    . $Bio::Parse::Format::FASTA3::RANK_START
    . "|"
    . $Bio::Parse::Format::FASTA3X::RANK_START
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::RANK_START
    . ")";
my $RANK_END = "(?:"
    . $Bio::Parse::Format::FASTA1::RANK_END
    . "|"
    . $Bio::Parse::Format::FASTA2::RANK_END
    . "|"
    . $Bio::Parse::Format::FASTA3::RANK_END
    . "|"
    . $Bio::Parse::Format::FASTA3X::RANK_END
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::RANK_END
    . ")";

my $MATCH_START = "(?:"
    . $Bio::Parse::Format::FASTA1::MATCH_START
    . "|"
    . $Bio::Parse::Format::FASTA2::MATCH_START
    . "|"
    . $Bio::Parse::Format::FASTA3::MATCH_START
    . "|"
    . $Bio::Parse::Format::FASTA3X::MATCH_START
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::MATCH_START
    . ")";
my $MATCH_END = "(?:"
    . $Bio::Parse::Format::FASTA1::MATCH_END
    . "|"
    . $Bio::Parse::Format::FASTA2::MATCH_END
    . "|"
    . $Bio::Parse::Format::FASTA3::MATCH_END
    . "|"
    . $Bio::Parse::Format::FASTA3X::MATCH_END
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::MATCH_END
    . ")";

my $SUM_START = "(?:"
    . $Bio::Parse::Format::FASTA1::SUM_START
    . "|"
    . $Bio::Parse::Format::FASTA2::SUM_START
    . "|"
    . $Bio::Parse::Format::FASTA3::SUM_START
    . "|"
    . $Bio::Parse::Format::FASTA3X::SUM_START
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::SUM_START
    . ")";
my $SUM_END = "(?:"
    . $Bio::Parse::Format::FASTA1::SUM_END
    . "|"
    . $Bio::Parse::Format::FASTA2::SUM_END
    . "|"
    . $Bio::Parse::Format::FASTA3::SUM_END
    . "|"
    . $Bio::Parse::Format::FASTA3X::SUM_END
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::SUM_END
    . ")";

my $ALN_START = "(?:"
    . $Bio::Parse::Format::FASTA1::ALN_START
    . "|"
    . $Bio::Parse::Format::FASTA2::ALN_START
    . "|"
    . $Bio::Parse::Format::FASTA3::ALN_START
    . "|"
    . $Bio::Parse::Format::FASTA3X::ALN_START
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::ALN_START
    . ")";
my $ALN_END = "(?:"
    . $Bio::Parse::Format::FASTA1::ALN_END
    . "|"
    . $Bio::Parse::Format::FASTA2::ALN_END
    . "|"
    . $Bio::Parse::Format::FASTA3::ALN_END
    . "|"
    . $Bio::Parse::Format::FASTA3X::ALN_END
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::ALN_END
    . ")";

my $TRAILER_START = "(?:"
    . $Bio::Parse::Format::FASTA1::TRAILER_START
    . "|"
    . $Bio::Parse::Format::FASTA2::TRAILER_START
    . "|"
    . $Bio::Parse::Format::FASTA3::TRAILER_START
    . "|"
    . $Bio::Parse::Format::FASTA3X::TRAILER_START
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::TRAILER_START
    . ")";
my $TRAILER_END = "(?:"
    . $Bio::Parse::Format::FASTA1::TRAILER_END
    . "|"
    . $Bio::Parse::Format::FASTA2::TRAILER_END
    . "|"
    . $Bio::Parse::Format::FASTA3::TRAILER_END
    . "|"
    . $Bio::Parse::Format::FASTA3X::TRAILER_END
    . "|"
    . $Bio::Parse::Format::FASTA2_GCG::TRAILER_END
    . ")";


#Generic get_entry() and new() constructors for all FASTA style parsers:
#determine program and version and coerce appropriate subclass.

my $DEBUG = 0;

#Consume one entry-worth of input on text stream associated with $file and
#return a new FASTA instance.
sub get_entry {
    my $text = shift;
    my $line = '';
    my $data = 0;

    my ($type, $prog, $version, $format) = ('Bio::Parse::Format::FASTA');
    my $GCG = 0;

    while ($text->getline(\$line)) {

        #start of entry
        if ($line =~ /$ENTRY_START/o and !$data) {
            $text->start_count();
            $data = 1;
            #fall through for version tests
        }

        #escape iteration if we've hit the alignment section
        next  if $line =~ /$ALN_START/;

        #try to determine program from header on stderr

        if ($line =~ /^\s*(\S+)\s+(?:searches|compares|produces|translates|performs)/) {
            $prog = $1;
            next;
        }

        #try to determine version from header on stderr

        if ($line =~ /^\s*version\s+(\d)\./) {
            $version = $1;   #eg, version 3.4
            next;
        }

        if ($line =~ /^\s*version\s+(3\d)/) {
            $version = '3X'; #eg, version 34
            next;
        }

        if ($line =~ /^\s*v(\d+)\.\d+\S\d+/) {
            $version = $1;
            next;
        }

        #try to determine program and version by minor differences on stdout

        if ($line =~ /The best scores are:\s+initn\s+init1\s+opt\s*$/) {
            $prog    = 'FASTA'    unless defined $prog;    #guess!
            $version = 1          unless defined $version;
            next;
        }

        if ($line =~ /The best scores are:\s+initn\s+init1\s+opt\s+z-sc/) {
            #matches FASTA2,FASTA3,TFASTX3, but next rules commit first
            $prog    = 'FASTA'    unless defined $prog;    #guess!
            $version = 2          unless defined $version;
            next;
        }

        if ($line =~ /The best scores are:\s+init1\s+initn\s+opt\s+z-sc/) {
            #matches GCG FASTA2
            $prog    = 'FASTA'    unless defined $prog;    #guess!
            $version = 2          unless defined $version;
            $GCG = 1;
            next;
        }

        if ($line =~ /ALIGN calculates/) {
            #matches FASTA2/ALIGN
            $prog    = 'ALIGN'     unless defined $prog;
            $version = 2           unless defined $version;
            next;
        }

        if ($line =~ /LALIGN finds/) {
            #matches FASTA2/LALIGN with stderr
            $prog    = 'LALIGN'    unless defined $prog;
            $version = 2           unless defined $version;
            next;
        }

        if ($line =~ /Comparison of:/) {
            #matches FASTA2/LALIGN
            $prog    = 'LALIGN'    unless defined $prog;
            $version = 2           unless defined $version;
            next;
        }

        if ($line =~ /SSEARCH searches/) {
            #matches FASTA2/SSEARCH
            $prog    = 'SSEARCH'   unless defined $prog;
            $version = 2           unless defined $version;
            next;
        }

        if ($line =~ /^(\S+)\s+\((\d+)/) {
            $prog    = $1          unless defined $prog;
            $version = $2          unless defined $version;
            next;
        }

        if ($line =~ /frame-shift:/) {
            #guess
            $prog    = 'TFASTX';
            next;
        }

        if ($line =~ /$GCG_JUNK/) {
            #matches GCG
            $GCG = 1;
            next;
        }

    }
    return 0  unless $data;

    unless (defined $prog and defined $version) {
        die "get_entry() top-level FASTA parser could not determine program/version\n";
    }

    unless (exists $VERSIONS{$version} and
            grep(/^$prog$/i, @{$VERSIONS{$version}}) > 0) {
        die "get_entry() parser for program '$prog' version '$version' not implemented\n";
    }

    $prog    = lc $prog;
    $format  = uc $prog;
    $version =~ s/-/_/g;

    if ($GCG) {
        $type = "Bio::Parse::Format::FASTA${version}_GCG::$prog";
    } else {
        $type = "Bio::Parse::Format::FASTA${version}::$prog";
    }

    #reuse packages: tfastx,tfasty,tfastxy -> tfasta
    $type =~ s/::tfast([xy]|xy)$/::tfasta/;

    #reuse packages: fasty -> fastx
    $type =~ s/::fasty$/::fastx/;

    #reuse packages: fastf,fasts -> fastm
    $type =~ s/::fast[fs]$/::fastm/;

    #reuse packages: tfastf,tfasts -> tfastm
    $type =~ s/::tfast[fs]$/::tfastm/;

    #warn "\nprog= $prog  version= $version (GCG=$GCG) type= $type\n";

    load_parser_class($type);

    my $self = $type->new(undef, $text, $text->get_start(), $text->get_stop());

    $self->{'format'}  = $format;
    $self->{'version'} = $version;

    $self;
}

sub load_parser_class {
    my $class = shift;
    $class =~ s/::/\//g;
    require "$class.pm";
}

#Parse one entry: generic for all FASTA[12]
#(FASTA3 $MATCH_START definition conflicts)
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
            $scan->scan_until($RANK_END);
            $self->push_record('RANK',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #Hit lines
        if ($line =~ /$MATCH_START/o) {
            $scan->scan_until($MATCH_END);
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

sub parse_frame {
    return 'f'  unless defined $_[0]; #default is forwards
    my $s = "$_[0]";
    return 'f'  if $s =~ /^\s*$/;     #default is forwards
    return 'f'  if $s eq 'f';
    return 'r'  if $s eq 'r';
    return 'r'  if $s eq 'rev-comp';
    return $s   if $s =~ /^[123]$/;   #seen in: tfasta_3.4t23
    return $s   if $s =~ /^[456]$/;   #seen in: tfasta_3.4t23
    return 'F'  if $s eq 'F';         #seen in: tfastx2.0u
    return 'R'  if $s eq 'R';         #seen in: tfastx2.0u
    return $s; #silently accept whatever was given
}

sub parse_orient {
    return '+'  unless defined $_[0]; #default is forwards
    my $s = "$_[0]";
    return '+'  if $s eq 'f';
    return '-'  if $s eq 'r';
    return '-'  if $s eq 'rev-comp';
    return '+'  if $s =~ /^[123]$/;   #seen in: tfasta_3.4t23
    return '-'  if $s =~ /^[456]$/;   #seen in: tfasta_3.4t23
    return '+'  if $s eq 'F';         #seen in: tfastx2.0u
    return '-'  if $s eq 'R';         #seen in: tfastx2.0u
    warn "WARNING: parse_orient: unrecognised value '$s'\n";
    return '?'; #what orientation?!!
}


###########################################################################
package Bio::Parse::Format::FASTA::HEADER;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new { die "$_[0]::new() virtual function called\n" }

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $field (sort $self->list_attrs) {
        $s .= sprintf "$x%20s -> %s\n", $field,  $self->{$field};
    }
    return $s;
}


###########################################################################
package Bio::Parse::Format::FASTA::RANK;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new { die "$_[0]::new() virtual function called\n" }

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $hit (@{$self->{'hit'}}) {
        foreach my $field (sort keys %$hit) {
            $s .= sprintf "$x%20s -> %s\n", $field, $hit->{$field};
        }
    }
    return $s;
}


###########################################################################
package Bio::Parse::Format::FASTA::TRAILER;

use vars qw(@ISA);

@ISA   = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);

    $self->{'trailer'} = $scan->read_remainder();

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $field (sort $self->list_attrs) {
        $s .= sprintf "$x%20s -> %s\n", $field,  $self->{$field};
    }
    return $s;
}


###########################################################################
package Bio::Parse::Format::FASTA::MATCH;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #identifier lines
        if ($line =~ /$SUM_START/o) {
            $scan->scan_until_inclusive($SUM_END);
            $self->push_record('SUM',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #fragment hits: terminated by several possibilities
        if ($line =~ /$ALN_START/o) {
            $scan->scan_until($ALN_END);
            $self->push_record('ALN',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #blank line or empty record: ignore
        next    if $line =~ /$NULL/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $field (sort $self->list_attrs) {
        $s .= sprintf "$x%20s -> %s\n", $field,  $self->{$field};
    }
    return $s;
}


###########################################################################
package Bio::Parse::Format::FASTA::MATCH::SUM;

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new { die "$_[0]::new() virtual function called\n" }

sub get_sibling {
    my ($self, $num) = @_;
    return $self->{'parent'}->{'blockkeeper'}->get_block('SUM', $num)->{'record'};
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    foreach my $field (sort $self->list_attrs) {
        $s .= sprintf "$x%20s -> %s\n", $field,  $self->{$field};
    }
    return $s;
}


###########################################################################
package Bio::Parse::Format::FASTA::MATCH::ALN;

use Bio::Util::Math qw(min max);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub get_header {
    $_[0]->{'parent'}->{'parent'}->{'blockkeeper'}->get_block('HEADER')->{'record'};
}

sub get_summary {
    $_[0]->{'parent'}->{'blockkeeper'}->get_block('SUM')->{'record'};
}

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    my $keys_and_depth = sub {
        my ($query, $sbjct) = @_;
        #determine depth of sequence labels at left: take the
        #shorter since either sequence may begin with a gap
        $$query =~ /^(\s*\S+\s*)/o; my $x = length $1;
        $$sbjct =~ /^(\s*\S+\s*)/o; my $y = length $1;
        my $depth = min($x, $y);
        #warn "depth: $depth\n";
        #recover sequence row names
        my $qkey = substr($$query, 0, $depth);
        my $skey = substr($$sbjct, 0, $depth);
        return ($qkey, $skey, $depth);
    };

    my $append_ruler = sub {
        my ($ruler, $piece, $depth) = @_;
        if ($$ruler eq '') {
            $$ruler .= $$piece;
            return $$ruler;
        }
        my ($num, $junk, $rlen) = (0, $depth, length($$ruler));
        my $s = substr($$piece, 0, $depth);
        if ($s =~ /(\d+)/o) {
            $num  = length $1;
            $junk -= $num;
        }
        $ruler = substr($$ruler, 0, $rlen-$num) . substr($$piece, $junk);
        #warn "ruler: $num/$junk/$depth/$rlen [$ruler]\n";
        $ruler;
    };

    my $summary_info = sub {
        my ($self) = @_;
        #look in the parent MATCH record's already parsed SUM record
        my $sum = $self->get_summary;
        #get any alignment ranges (but these only cover the aligned
        #regions; terminal overhanging sequences are not counted!
        my ($q1, $q2, $s1, $s2) = (0, 0, 0, 0);
        if (exists $sum->{'ranges'} and
            $sum->{'ranges'} =~ /(\d+)-(\d+):(\d+)-(\d+)/) {
            ($q1, $q2, $s1, $s2) = ($1, $2, $3, $4);
        }
       #warn "\n>>> $sum->{'id'}\n"  if $DEBUG;
       #warn "ranges: [$q1, $q2], [$s1, $s2]\n\n"  if $DEBUG;
        [ [$q1, $q2], [$s1, $s2] ];
    };

    my $start_info = sub {
        my ($self, $name, $ruler, $string, $depth, $rstart) = @_;
        my ($ai, $si, $ni, $num, $usedr) = (0, 0, 0, 0, 0);

        #find the index of the start of sequence
        if ($$string =~ /(^[-\s]*)\S/o) {
            $si = length($1);
        } else {
            $self->die("no sequence:", $name);
        }

        #find the index of the first ruler number not over a gap (it happens)
        my $rcopy = $$ruler;
        my $pos = -$depth - 1;
        while ($rcopy ne '') {
            if ($rcopy =~ /(\s*(\d+))/o) {
                $pos += length $1;
                my $tmp = $2;
                my $ch = substr($$string, $pos, 1);
                #warn "look: [$rcopy]  $pos  [$ch]  $tmp\n";
                $rcopy = $'; #'emacs;
                if ($ch !~ /[-\/\\]/o) {
                    ($ni, $num) = ($pos, $tmp);
                    #warn "exit: [$rcopy]  $pos  [$ch]  $tmp\n";
                    last;
                }
                next;
            }
            last;
        }
        #give up and take any range numbers
        if ($ni == 0 and $num == 0) {
            ($ni, $num, $usedr) = ($si, $rstart, 1);
        }

        #count gaps in the fragment between seqstart and label
        my $gc = 0;
        if ($si < $ni) {  #seqstart ... label
            my $s = substr($$string, $si, $ni-$si);
            $gc = $s =~ tr/-\/\\//;
        }
        [$rstart, $ai, $si, $ni, $num, $gc, $usedr];
    };

    my $stop_info = sub {
        my ($self, $name, $ruler, $string, $depth, $rend) = @_;
        my ($ai, $si, $ni, $num, $usedr) = (0, 0, 0, 0, 0);

        #find the index of the end of sequence
        if ($$string =~ /\S([-\s]*)$/o) {
            $ai = length($$string) -1;
            $si = $ai - length($1);
        } else {
            $self->die("no sequence:", $name);
        }

        #find the index of the last ruler number not over a gap (it happens)
        my $rcopy = $$ruler;
        my $pos = length($rcopy) - $depth - 1;
        while ($rcopy ne '') {
            if ($rcopy =~ /(\d+)(\s*)$/o) {
                $pos -= length $2;
                my $tmp = $1;
                my $ch = substr($$string, $pos, 1);
                #warn "look: [$rcopy]  $pos  [$ch]  $tmp\n";
                $rcopy = $`;
                if ($ch !~ /[-\/\\]/o) {
                    ($ni, $num) = ($pos, $tmp);
                    #warn "exit: [$rcopy]  $pos  [$ch]  $tmp\n";
                    last;
                }
                $pos -= length $tmp;
                next;
            }
            last;
        }
        #give up and take any range numbers
        if ($ni == 0 and $num == 0) {
            ($ni, $num, $usedr) = ($si, $rend, 1);
        }

        #count non-sequence in the fragment between label and seqend
        my $gc = 0;
        if ($ni < $si) {  #label ... seqend
            my $s = substr($$string, $ni, $si-$ni);
            $gc = $s =~ tr/-\/\\//;
        }
        [$rend, $ai, $si, $ni, $num, $gc, $usedr];
    };

    my $align_info = sub {
        my ($self, $sum, $name, $ruler, $string, $depth) = @_;
        my $r1  = $sum->[0];  #summary from
        my $r2  = $sum->[1];  #summary to
        my $start = &$start_info($self, $name, $ruler, $string, $depth, $r1);
        my $stop  = &$stop_info($self, $name, $ruler, $string, $depth, $r2);
        #count frameshifts lying between the extreme ruler numbers
        #my ($n1, $n2) = ($start->[2], $stop->[2]);
        #my $s = substr($$string, $n1, $n2-$n1+1);
        #count frameshifts of types '/' and '\'
        my $s = $$string;
        my $fsf = $s =~ tr/[\/]//;
        my $fsb = $s =~ tr/[\\]//;
        # if ($DEBUG) {
        #     warn "[$$ruler]\n[@{['.'x$depth]}$$string]\n";
        #     warn sprintf("start:  [r1:%4d  a:%4d  b:%4d  x:%4d  X:%4d  gc: %d  usedr: %d]\n", @$start);
        #     warn sprintf("stop:   [r2:%4d  a:%4d  b:%4d  x:%4d  Y:%4d  gc: %d  usedr: %d]\n", @$stop);
        #     warn "fshift: [/ $fsf, \\ $fsb]\n\n";
        # }
        [ $start, $stop, $fsf, $fsb ];
    };

    my ($qrule, $query, $align, $sbjct, $srule) = ('', '', '', '', '');
    my ($first_pass, $depth, $qkey, $skey) = (1, 0, '', '');

    #read records
    while (defined ($line = $scan->read_line(1))) {
        my @tmp = ();

        #initial empty line before the alignment ruler
        next if $line =~ /^$/;

        #warn "$first_pass, [$line]\n";

        if ($line =~ /^\s+(?:\d+)?/) {
            #warn "QUERY RULER\n";
            $tmp[0] = $line;                #query ruler
            $tmp[1] = $scan->read_line(1);  #query sequence
            $tmp[2] = $scan->read_line(1);  #match pattern
            $tmp[3] = $scan->read_line(1);  #sbjct sequence
            $tmp[4] = $scan->read_line(1);  #sbjct ruler

            if ($first_pass) {
                ($qkey, $skey, $depth) = &$keys_and_depth(\$tmp[1], \$tmp[3]);
                #warn "ENTRY [$qkey], [$skey]\n";
                $first_pass = 0;
            }

        } elsif (index($line, $qkey) == 0) {
            #warn "QUERY (####)\n";
            $tmp[1] = $line;                #query sequence
            $tmp[2] = $scan->read_line(1);  #match pattern
            $tmp[3] = $scan->read_line(1);  #sbjct sequence
            $tmp[4] = $scan->read_line(1);  #sbjct ruler

        } elsif (index($line, $skey) == 0) {
            #warn "SBJCT (####)\n";
            $tmp[3] = $line;                #sbjct sequence
            $tmp[4] = $scan->read_line(1);  #sbjct ruler

        } else {
            $self->die("unexpected line: [$line]\n");
        }
        map { $tmp[$_] = '' unless defined $tmp[$_] } 0..@tmp-1;

        #warn "#0# [$tmp[0]]\n";
        #warn "#1# [$tmp[1]]\n";
        #warn "#2# [$tmp[2]]\n";
        #warn "#3# [$tmp[3]]\n";
        #warn "#4# [$tmp[4]]\n";

        #pad query/match/sbjct lines
        my $len = max(length($tmp[1]), length($tmp[3]));
        $tmp[0] .= ' ' x ($len-length $tmp[0]);
        $tmp[1] .= ' ' x ($len-length $tmp[1]);
        $tmp[2] .= ' ' x ($len-length $tmp[2]);
        $tmp[3] .= ' ' x ($len-length $tmp[3]);
        $tmp[4] .= ' ' x ($len-length $tmp[4]);

        #strip leading name from alignment rows
        $tmp[1] = substr($tmp[1], $depth);
        $tmp[2] = substr($tmp[2], $depth);
        $tmp[3] = substr($tmp[3], $depth);

        #grow the ruler
        $qrule = &$append_ruler(\$qrule, \$tmp[0], $depth);
        $srule = &$append_ruler(\$srule, \$tmp[4], $depth);

        #grow the alignment
        $query .= $tmp[1];
        $align .= $tmp[2];
        $sbjct .= $tmp[3];
    }

    #query/sbjct start/stop positions from parent summary
    my $sum = &$summary_info($self);

    #query/sbjct start/stop positions and numbering from ruler
    my $qinfo = &$align_info($self, $sum->[0], $qkey, \$qrule, \$query, $depth);
    my $sinfo = &$align_info($self, $sum->[1], $skey, \$srule, \$sbjct, $depth);

    #query/sbjct orientations
    my ($query_orient, $qorient) = $self->get_query_orient($qinfo);
    my ($sbjct_orient, $sorient) = $self->get_sbjct_orient($sinfo);

    #query/sbjct recomputed start/stop
    my ($query_start, $query_stop) =
        $self->get_start_stop('qry', $qorient, $qinfo, $self->query_base);
    my ($sbjct_start, $sbjct_stop) =
        $self->get_start_stop('hit', $sorient, $sinfo, $self->sbjct_base);

   #warn "Final: [$query_orient $query_start,$query_stop] [$sbjct_orient $sbjct_start $sbjct_stop]\n\n"  if $DEBUG;

    my $x;

    #query length
    $x = $query;
    $x =~ tr/- //d;
    my $query_length = length($x);

    #sbjct length
    $x = $sbjct;
    $x =~ tr/- //d;
    my $sbjct_length = length($x);

    #query_leader
    $query =~ /^(\s*)/;
    my $query_leader = length $1;

    #query_trailer
    $query =~ /(\s*)$/;
    my $query_trailer = length $1;

    #sbjct_leader
    $sbjct =~ /^(\s*)/;
    my $sbjct_leader = length $1;

    #sbjct_trailer
    $sbjct =~ /(\s*)$/;
    my $sbjct_trailer = length $1;

    $self->{'query'} = $query;
    $self->{'align'} = $align;
    $self->{'sbjct'} = $sbjct;

    #warn "QUERY ($query_start, $query_stop, $query_length, $query_leader)\n";
    #warn "SBJCT ($sbjct_start, $sbjct_stop, $sbjct_length, $sbjct_leader)\n";
    #warn "EXIT  ($query_orient, $query_start, $query_stop) ($sbjct_orient, $sbjct_start, $sbjct_stop)\n";

    $self->{'query_orient'}  = $query_orient;
    $self->{'query_start'}   = $query_start;
    $self->{'query_stop'}    = $query_stop;
    $self->{'query_leader'}  = $query_leader;
    $self->{'query_trailer'} = $query_trailer;

    $self->{'sbjct_orient'}  = $sbjct_orient;
    $self->{'sbjct_start'}   = $sbjct_start;
    $self->{'sbjct_stop'}    = $sbjct_stop;
    $self->{'sbjct_leader'}  = $sbjct_leader;
    $self->{'sbjct_trailer'} = $sbjct_trailer;

    $self;
}

#summary orientation may refer to query or to subjct and may differ from
#sequence numbering, depending on program and version. so get both.
sub get_query_orient {
    my ($self, $info) = @_;
    #ruler orientation
    my $r1 = $info->[0][0]; my $r2 = $info->[1][0];
    my $n1 = $info->[0][4]; my $n2 = $info->[1][4];
    #label orientation: child can override method
    my $lo = $self->query_orient;
    my $ro;
    if ($r1 != $r2) {
        $ro = $r2 < $r1 ? '-' : '+';
    } elsif ($n1 != $n2) {
        $ro = $n2 < $n1 ? '-' : '+';
    } else {
        $ro = $lo;
    }
    if ($lo ne $ro) {
        my $id = $self->get_summary->{'id'};
        warn "WARNING: query orientation conflict for '$id': summary says '$lo', but range is '$ro'\n"  if $DEBUG;
    }
    ($lo, $ro);
}

sub get_sbjct_orient {
    my ($self, $info) = @_;
    #ruler orientation
    my $r1 = $info->[0][0]; my $r2 = $info->[1][0];
    my $n1 = $info->[0][4]; my $n2 = $info->[1][4];
    #label orientation: child can override method
    my $lo = $self->sbjct_orient;
    my $ro;
    if ($r1 != $r2) {
        $ro = $r2 < $r1 ? '-' : '+';
    } elsif ($n1 != $n2) {
        $ro = $n2 < $n1 ? '-' : '+';
    } else {
        $ro = $lo;
    }
    if ($lo ne $ro) {
        my $id = $self->get_summary->{'id'};
        warn "WARNING: sbjct orientation conflict for '$id': summary says '$lo', but range is '$ro'\n"  if $DEBUG;
    }
    ($lo, $ro);
}

##Generic functions to determine orientations: Override as needed.
sub query_orient { '+' }
sub sbjct_orient { '+' }

#Generic functions to fix start/stop positions. Override as needed.
sub query_base { return 1 }
sub sbjct_base { return 1 }

#subclasses override this if they need to do less work
sub get_start_stop {
    my ($self, $tgt, $orient, $info, $base) = @_;

    my ($start_info, $stop_info, $fsf, $fsb) = @$info;

   #warn ">>>>>> $tgt, $orient\n"  if $DEBUG;

    #              7            20
    #  ---------QWERTY  //  UIOPASDF----------
    #  a        b  x    //       y c         d

    # a:     alignment start: should be 0
    # b:     sequence  start: >= 0
    # c:     sequence  stop:  >= b
    # d:     alignment stop:  should be length(sequence)-1
    # x:     numbering start: >= 0  (could be < b?)
    # y:     numbering stop:  >= c
    # nx:    number start
    # ny:    number stop
    # gc1/2: count of non-sequence characters before/after the number
    # r1/2:  summary range start/stop number

    my ($r1, $a, $b, $x, $nx, $gc1, $usedr1) = @{$start_info};
    my ($r2, $d, $c, $y, $ny, $gc2, $usedr2) = @{$stop_info};

    my $delta1 = abs($x - $b) - $gc1;
    my $delta2 = abs($y - $c) - $gc2;

    if ($x < $b) {
        $delta1 = -$delta1;
    }
    if ($y > $c) {
        $delta2 = -$delta2;
    }

    my ($start, $stop, $shift) = (0, 0, 0);

    #extrapolate endpoints from the ruler numbers by the above deltas:
    #- scale by $base (1=protein or 3=DNA)
    #- account for the rest of the codon at the end (+/- 2) if in base 3
    #- shift according to the alignment summary range values
   #warn "$tgt(i): $orient,$base bc= $b,$c  xy= $x,$y  nxy= $nx,$ny  delta1/2= $delta1,$delta2  usedr= $usedr2\n"  if $DEBUG;
    if ($orient eq '+') {
        if ($base == 1) {
            $start = $nx - $delta1;
            $stop  = $ny + $delta2;
        } else {
            $start = $nx - $base * $delta1;
            $stop  = $ny + $base * $delta2;
            $stop += 2  unless $usedr2;
            if ($r1 != 0) {
                $shift = $r1 - $start;
                $start -= $shift;
                $stop  -= $shift;
            }
           #warn "$tgt(i): $orient,$base start/stop: $start,$stop  shift: $shift\n"  if $DEBUG;
        }
    } else {
        if ($base == 1) {
            $start = $nx + $delta1;
            $stop  = $ny - $delta2;
        } else {
            $start = $nx + $base * $delta1;
            $stop  = $ny - $base * $delta2;
            $stop -= 2  unless $usedr2;
            if ($r1 != 0) {
                $shift = $r1 - $start;
                $start += $shift;
                $stop  += $shift;
            }
           #warn "$tgt(i): $orient,$base start/stop: $start,$stop  shift: $shift\n"  if $DEBUG;
        }
    }
   #warn "$tgt(o): $orient,$base start/stop: $start,$stop\n"  if $DEBUG;
   #warn "\n"  if $DEBUG;

    return ($start, $stop);
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> '%s'\n", 'query',          $self->{'query'};
    $s .= sprintf "$x%20s -> '%s'\n", 'align',          $self->{'align'};
    $s .= sprintf "$x%20s -> '%s'\n", 'sbjct',          $self->{'sbjct'};
    $s .= sprintf "$x%20s -> %s\n",   'query_orient',   $self->{'query_orient'};
    $s .= sprintf "$x%20s -> %s\n",   'query_start',    $self->{'query_start'};
    $s .= sprintf "$x%20s -> %s\n",   'query_stop',     $self->{'query_stop'};
    $s .= sprintf "$x%20s -> %s\n",   'query_leader',   $self->{'query_leader'};
    $s .= sprintf "$x%20s -> %s\n",   'query_trailer',  $self->{'query_trailer'};
    $s .= sprintf "$x%20s -> %s\n",   'sbjct_orient',   $self->{'sbjct_orient'};
    $s .= sprintf "$x%20s -> %s\n",   'sbjct_start',    $self->{'sbjct_start'};
    $s .= sprintf "$x%20s -> %s\n",   'sbjct_stop' ,    $self->{'sbjct_stop'};
    $s .= sprintf "$x%20s -> %s\n",   'sbjct_leader',   $self->{'sbjct_leader'};
    $s .= sprintf "$x%20s -> %s\n",   'sbjct_trailer',  $self->{'sbjct_trailer'};
    return $s;
}


###########################################################################
1;
