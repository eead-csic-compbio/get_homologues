# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# Handles: fasta  1.x
#
# Acknowledgements.
#
# Christophe Leroy, 22/8/97 for original fasta1 parser.
#
###########################################################################
package Bio::Parse::Format::FASTA1;

use Bio::Parse::Format::FASTA;
use strict;

use vars qw(@ISA

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
             '1' => [
                     'FASTA',
                    ],
            );

$NULL  = '^\s*$';

$ENTRY_START   = '(?:'
    . '^\s*fasta.*searches a sequence data bank'
    . '|'
    . '^\s*\S+\s*[,:]\s+\d+\s+(?:aa|nt)'
    . ')';
$ENTRY_END     = 'Library scan:';

$HEADER_START  = $ENTRY_START;
$HEADER_END    = '^The best scores are:';

$RANK_START    = $HEADER_END;
$RANK_END      = $NULL;

$TRAILER_START = $ENTRY_END;
$TRAILER_END   = $ENTRY_END;

$MATCH_START   = '^\S{7}.*\d+\s+\d+\s+\d+\s*$';
$MATCH_END     = "(?:$MATCH_START|$ENTRY_END)";

$SUM_START     = $MATCH_START;
$SUM_END       = $NULL;

$ALN_START     = '^(?:\s+\d+\s+|\s+$)';  #the ruler
$ALN_END       = $MATCH_END;

#Parse one entry
sub new { my $self=shift; $self->SUPER::new(@_) }


###########################################################################
package Bio::Parse::Format::FASTA1::MATCH::ALN;

use Bio::Util::Math qw(max);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Format::FASTA::MATCH::ALN);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    my ($query_start, $query_stop, $sbjct_start, $sbjct_stop) = (0,0,0,0);
    my ($query_leader,$query_trailer,$sbjct_leader,$sbjct_trailer) = (0,0,0,0);
    my ($query_length, $sbjct_length, $query_orient, $sbjct_orient);
    my ($query, $align, $sbjct, $len) = ('','','');
    my ($x, $querykey, $sbjctkey, $depth);

    #first record
    my @tmp = ();

    #warn $scan->inspect_stream;

    $line = $scan->read_line(1);

    #initial empty line before the alignment ruler
    if ($line =~ /^$/) {
        $line = $scan->read_line(1);
    }

    if ($line =~ /^\s*(\S+)\s+([^0-9 ]+)$/) {
        #no query ruler on first line, due to pure leading gaps seen
        #in ALIGN and ALIGN0 or very short sequences in LALIGN.
        #very short sequences can mean that start/stop is missing!
        $tmp[0] = '';                     #missing ruler
        $tmp[1] = $line;                  #query sequence
        $tmp[2] = $scan->read_line(1);    #match pattern
        $tmp[3] = $scan->read_line(1);    #sbjct sequence
        $tmp[4] = $scan->read_line(1);    #sbjct ruler
    } elsif ($line =~ /^\s*\d+/) {
        #query ruler present on first line
        $tmp[0] = $line;                  #query ruler
        $tmp[1] = $scan->read_line(1);    #query sequence
        $tmp[2] = $scan->read_line(1);    #match pattern
        $tmp[3] = $scan->read_line(1);    #sbjct sequence
        $tmp[4] = $scan->read_line(1);    #sbjct ruler
    } elsif ($line =~ /^\s+$/) {
        #ruler has no numbers as sequence is too short
        $tmp[0] = '';                     #missing ruler
        $tmp[1] = $scan->read_line(1);    #query sequence
        $tmp[2] = $scan->read_line(1);    #match pattern
        $tmp[3] = $scan->read_line(1);    #sbjct sequence
        $tmp[4] = $scan->read_line(1);    #sbjct ruler
    } else {
        $self->die("unexpected line: [$line]\n");
    }

    #determine depth of sequence labels at left: take the shorter
    #since either sequence may begin with a gap
    $tmp[1] =~ /^(\s*\S+\s*)/;
    $depth = length $1;
    $tmp[3] =~ /^(\s*\S+\s*)/;
    $x = length $1;
    $depth = ($depth < $x ? $depth : $x);

    #recover sequence row names
    $querykey = substr($tmp[1], 0, $depth);
    $sbjctkey = substr($tmp[3], 0, $depth);

    #strip leading name
    $tmp[1] = substr($tmp[1], $depth);
    $tmp[2] = substr($tmp[2], $depth);
    $tmp[3] = substr($tmp[3], $depth);

    #warn "#0##$tmp[0]\n";
    #warn "#1##$tmp[1]\n";
    #warn "#2##$tmp[2]\n";
    #warn "#3##$tmp[3]\n";
    #warn "#4##$tmp[4]\n";

    #query/sbjct/match lines
    #warn "|$tmp[1]|\n|$tmp[2]|\n|$tmp[3]|\n";
    $len = max(length($tmp[1]), length($tmp[3]));
    $tmp[1] .= ' ' x ($len-length $tmp[1])    if length $tmp[1] < $len;
    $tmp[2] .= ' ' x ($len-length $tmp[2])    if length $tmp[2] < $len;
    $tmp[3] .= ' ' x ($len-length $tmp[3])    if length $tmp[3] < $len;
    $query = $tmp[1];
    $align = $tmp[2];
    $sbjct = $tmp[3];

    #initialise query length
    $x = $tmp[1]; $x =~ tr/- //d; $query_length = length($x);

    #first query start/stop, if any
    if ($tmp[0] =~ /^\s*(\d+)/) {
        $query_start = $1;
    }
    if ($tmp[0] =~ /(\d+)\s*$/) {
        $query_stop  = $1;
    }

    #compute length of query before query_start label
    $x = index($tmp[0], $query_start) + length($query_start) -1;
    $x = substr($tmp[1], 0, $x-$depth);  #whole leader
    $x =~ tr/- //d;                      #leader minus gaps
    $query_leader = length($x);

    #initialise sbjct length
    $x = $tmp[3]; $x =~ tr/- //d; $sbjct_length = length($x);

    #first sbjct start/stop, if any
    if ($tmp[4] =~ /^\s*(\d+)/) {
        $sbjct_start = $1;
    }
    if ($tmp[4] =~ /(\d+)\s*$/) {
        $sbjct_stop  = $1;
    }

    #compute length of sbjct before sbjct_start label
    $x = index($tmp[4], $sbjct_start) + length($sbjct_start) -1;
    $x = substr($tmp[3], 0, $x-$depth);  #whole leader
    $x =~ tr/- //d;                      #leader minus gaps
    $sbjct_leader = length($x);

    #warn "ENTRY ($querykey, $sbjctkey)\n";
    #warn "($query_length/$query_leader) ($sbjct_length/$sbjct_leader)\n";

    #remaining records
    while (defined ($line = $scan->read_line(1))) {

        next if $line =~ /^\s*$/;

        @tmp = ();

        if ($line =~ /^\s*\d+(?:\s|$)/) {
            #query ruler
            #warn "QUERY+RULER\n";

            $tmp[0] = $line;                     #query ruler

            $line = $scan->read_line(1);
            $tmp[1] = substr($line, $depth);     #query sequence

            $line = $scan->read_line(1);
            if ($line) {
                $tmp[2] = substr($line, $depth); #match pattern
            } else {
                $tmp[2] = ' ' x length $tmp[1];
            }

            $line = $scan->read_line(1);
            if ($line) {
                $tmp[3] = substr($line, $depth); #sbjct sequence
            } else {
                $tmp[3] = ' ' x length $tmp[1];
            }

            $line = $scan->read_line(1);
            if ($line) {
                $tmp[4] = $line;                 #sbjct ruler
            } else {
                $tmp[4] = '';
            }

        } elsif (index($line, $sbjctkey) == 0) {
            #sbjct sequence
            #warn "SBJCT (####)\n";

            if ($line) {
                $tmp[3] = substr($line, $depth); #sbjct sequence
            } else {
                $tmp[3] = ' ' x length $tmp[1];
            }

            $line = $scan->read_line(1);
            if ($line) {
                $tmp[4] = $line;                 #sbjct ruler
            } else {
                $tmp[4] = '';
            }

            $tmp[0] = '';
            $tmp[1] = ' ' x length $tmp[3];
            $tmp[2] = ' ' x length $tmp[3];

        } elsif (index($line, $querykey) == 0) {
            #query sequence (no ruler)
            #warn "QUERY (####)\n";

            $tmp[1] = substr($line, $depth);     #query sequence

            $line = $scan->read_line(1);
            if ($line) {
                $tmp[2] = substr($line, $depth); #match pattern
            } else {
                $tmp[2] = ' ' x length $tmp[1];
            }

            $line = $scan->read_line(1);
            if ($line) {
                $tmp[3] = substr($line, $depth); #sbjct sequence
            } else {
                $tmp[3] = ' ' x length $tmp[1];
            }

            $line = $scan->read_line(1);
            if ($line) {
                $tmp[4] = $line;                 #sbjct ruler
            } else {
                $tmp[4] = '';
            }

            $tmp[0] = '';

        } else {
            $self->die("unexpected line: [$line]\n");
        }

        #warn "#0##$tmp[0]\n";
        #warn "#1##$tmp[1]\n";
        #warn "#2##$tmp[2]\n";
        #warn "#3##$tmp[3]\n";
        #warn "#4##$tmp[4]\n";

        #subsequent query start/stop, if any
        if ($query_start < 1 and $tmp[0] =~ /^\s*(\d+)/) {
            $query_start = $1;
            #compute length of query before query_start label
            $x = index($tmp[0], $query_start) + length($query_start) -1;
            $x = substr($tmp[1], 0, $x-$depth);  #whole leader
            $x =~ tr/- //d;                      #leader minus gaps
            $query_leader = length($x);
        }
        if ($tmp[0] =~ /(\d+)\s*$/) {
            $query_stop  = $1;    #always update stop
        }

        #subsequent sbjct start/stop, if any
        if ($sbjct_start < 1 and $tmp[4] =~ /^\s*(\d+)/) {
            $sbjct_start = $1;
        }
        if ($tmp[4] =~ /(\d+)\s*$/) {
            $sbjct_stop  = $1;    #always update stop
        }

        #increment query length
        $x = $tmp[1]; $x =~ tr/- //d; $query_length += length($x);

        #increment sbjct length
        $x = $tmp[3]; $x =~ tr/- //d; $sbjct_length += length($x);

        #query/sbjct/match lines
        #warn "|$tmp[1]|\n|$tmp[2]|\n|$tmp[3]|\n";
        $len = max(length($tmp[1]), length($tmp[3]));
        if (length $tmp[1] < $len) {
            $tmp[1] .= ' ' x ($len-length $tmp[1]);
        }
        if (length $tmp[2] < $len) {
            $tmp[2] .= ' ' x ($len-length $tmp[2]);
        }
        if (length $tmp[3] < $len) {
            $tmp[3] .= ' ' x ($len-length $tmp[3]);
        }
        $query .= $tmp[1];
        $align .= $tmp[2];
        $sbjct .= $tmp[3];

        #warn "($query_length) ($sbjct_length)\n";
    }

    #warn "$query\n$sbjct\n";
    #warn $scan->inspect_stream;

    #determine query orientation and start/stop
    #warn "QUERY ($query_start, $query_stop, $query_length, $query_leader)\n";
    ($query_orient, $query_start, $query_stop) =
        $self->adjust_counts($query_start, $query_stop, $query_length,
                 $query_leader, $self->query_base());

    #determine sbjct orientation and start/stop
    #warn "SBJCT ($sbjct_start, $sbjct_stop, $sbjct_length, $sbjct_leader)\n";
    ($sbjct_orient, $sbjct_start, $sbjct_stop) =
        $self->adjust_counts($sbjct_start, $sbjct_stop, $sbjct_length,
                             $sbjct_leader, $self->sbjct_base());

    #warn "EXIT ($query_orient, $query_start, $query_stop) ($sbjct_orient, $sbjct_start, $sbjct_stop)\n";

    #query_leader/query_trailer
    $query =~ /^(\s*)/;
    $query_leader   = length $1;
    $query =~ /(\s*)$/;
    $query_trailer  = length $1;

    #sbjct_leader/sbjct_trailer
    $sbjct =~ /^(\s*)/;
    $sbjct_leader   = length $1;
    $sbjct =~ /(\s*)$/;
    $sbjct_trailer  = length $1;

    $self->{'query'} = $query;
    $self->{'align'} = $align;
    $self->{'sbjct'} = $sbjct;

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

#Generic functions to fix start/stop positions suitable for any fasta member
#that uses untranslated sequences where the stated sequence ranges use the
#same numbering system as the displayed sequences. Must be overridden as
#needed:
#
# fasta       dna x dna  or  pro x pro  (ok, base=1)
# fasta[xy]   dna x pro                 (override query, base=3)
# tfast[axy]  pro x dna                 (override sbjct, base=3)
#
sub query_base { return 1 }
sub sbjct_base { return 1 }

#determine orientation and adjust start/stop
sub adjust_counts {
    my ($self, $start, $stop, $length, $leader, $base) = @_;
    my $orient = '+';

    $orient = '-' if $start > $stop;

    $leader *= $base;
    $length *= $base;

    #counting untranslated sequence
    if ($orient eq '+') {
        $start -= $leader;
        $stop   = $start + $length - 1;
    } else {
        $start += $leader;
        $stop   = $start - $length + 1;
    }

    return ($orient, $start, $stop);
}


###########################################################################
1;
