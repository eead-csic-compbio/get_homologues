# Copyright (C) 1997-2015 Nigel P. Brown
# $Id: Convert.pm $

######################################################################
package Bio::MView::Convert;

use strict;

my $DEF_PAD = '-';  #default terminal gap character
my $DEF_GAP = '-';  #default internal gap character

my %Known_Convert_Mode =
    (
     #name
     'plain'      => 1,
     'fasta'      => 1,
     'clustal'    => 1,
     'pir'        => 1,
     'msf'        => 1,
     'rdb'        => 1,
    );

sub check_convert_mode {
    if (defined $_[0]) {
	if (exists $Known_Convert_Mode{$_[0]}) {
	    return lc $_[0];
	}
    }
    return map { lc $_ } sort keys %Known_Convert_Mode;
}

sub new {
    my $type = shift;
    my $self = {};

    $self->{'build'}   = shift;
    $self->{'align'}   = shift;
    $self->{'moltype'} = shift;
    $self->{'pad'}     = defined $_[0] ? shift : $DEF_PAD;
    $self->{'gap'}     = defined $_[0] ? shift : $DEF_GAP;

    bless $self, $type;
    $self;
}

###########################################################################
my $PLAIN_ID_WIDTH = 20;  #default width for id' field

#return alignment in 'plain' format
sub plain {
    my ($self, $idw) = (@_, $PLAIN_ID_WIDTH);
    my ($bld, $aln, $s) = ($self->{'build'}, $self->{'align'}, '');
    foreach my $rid ($aln->visible_ids) {
        my $w = length($self->{'build'}->uid2row($rid)->cid);
        $idw = $w  if $w > $idw;
    }
    foreach my $rid ($aln->visible_ids) {
        #warn "[$rid]\n";
        $s .= $self->plain_row($rid, $idw);
    }
    \$s;
}

sub plain_row {
    my ($self, $rid, $idw) = @_;
    my $row = $self->{'build'}->uid2row($rid);
    my $title = sprintf("%-${idw}s", substr($row->cid, 0, $idw));
    my $sequence = $row->seq($self->{'pad'}, $self->{'gap'});
    $title . " " . $sequence . "\n";
}


###########################################################################
my $PEARSON_SEQ_WIDTH = 70;

#return alignment in Pearson/FASTA format
sub pearson {
    my $self = shift;
    my ($bld, $aln, $s) = ($self->{'build'}, $self->{'align'}, '');
    foreach my $rid ($aln->visible_ids) {
        #warn "[$rid]\n";
        $s .= $self->pearson_row($rid);
    }
    \$s;
}

sub pearson_row {
    my ($self, $rid) = @_;
    my $row = $self->{'build'}->uid2row($rid);

    my $head = sub {
	my $self = shift;
	my $s = ">";
	#my $d = $row->num;
	#$s .= ((defined $d and $d ne '') ? "$d;" : "query;");
	$s .= $row->cid;
    };

    my $desc = sub {
	my $row = shift;
	my ($s, $d) = ('');
	$d = $row->desc;  $s .= " $d"  if $d ne '';
	$d = $row->data;  $s .= " $d"  if $d ne '';
	$d = $row->posn1; $s .= " $d"  if $d ne '';
	$d = $row->posn2; $s .= " $d"  if $d ne '';
	$s . "\n";
    };

    my $sequence = sub {
	my ($row, $pad, $gap) = @_;
	my $seq = $row->seq($pad, $gap);
	my $len = length($seq);
	my $s = '';
	for (my $i=0; $i<$len; $i+=$PEARSON_SEQ_WIDTH) {
	    $s .= substr($seq, $i, $PEARSON_SEQ_WIDTH) . "\n";
	}
	$s;
    };

    &$head($row)
        . &$desc($row)
        . &$sequence($row, $self->{'pad'}, $self->{'gap'});
}

###########################################################################
my $PIR_SEQ_WIDTH = 60;

#return alignment in PIR format
sub pir {
    my $self = shift;
    my ($bld, $aln, $s) = ($self->{'build'}, $self->{'align'}, '');
    foreach my $rid ($aln->visible_ids) {
        #warn "[$rid]\n";
        $s .= $self->pir_row($rid);
    }
    \$s;
}

sub pir_row {
    my ($self, $rid) = @_;
    my $row = $self->{'build'}->uid2row($rid);

    my $head = sub {
	my ($row, $moltype) = (@_, 'aa');
	my $s = $moltype eq 'aa' ? ">P1;" : ">XX;";
	#my $d = $row->num;
	#$s .= ((defined $d and $d ne '') ? "$d;" : "query;");
	$s .= $row->cid . "\n";
    };

    my $desc = sub {
	my $row = shift;
	my ($s, $d) = ('');
	$d = $row->desc;  $s .= ($s eq '' ? $d : " $d")  if $d ne '';
	$d = $row->data;  $s .= ($s eq '' ? $d : " $d")  if $d ne '';
	$d = $row->posn1; $s .= ($s eq '' ? $d : " $d")  if $d ne '';
	$d = $row->posn2; $s .= ($s eq '' ? $d : " $d")  if $d ne '';
	$s = '.'  if $s eq '';
	$s;
    };

    my $sequence = sub {
	my ($row, $pad, $gap) = @_;
	my $seq = $row->seq($pad, $gap);
	my $len = length($seq);
	my $s = '';
	for (my $i=0; $i<$len; $i+=$PIR_SEQ_WIDTH) {
	    $s .= "\n" . substr($seq, $i, $PIR_SEQ_WIDTH);
	}
	$s .= "\n"    if length($s) % ($PIR_SEQ_WIDTH+1) < 1 and $s ne '';
	$s .= "*\n\n";
    };

    &$head($row, $self->{'moltype'})
        . &$desc($row)
        . &$sequence($row, $self->{'pad'}, $self->{'gap'});
}

###########################################################################
#return alignment in RDB table format
sub rdb {
    my $self = shift;
    my ($bld, $aln, $s) = ($self->{'build'}, $self->{'align'}, '');
    $s .= $bld->index2row(0)->rdb_row('attr') . "\n";
    $s .= $bld->index2row(0)->rdb_row('form') . "\n";
    foreach my $rid ($aln->visible_ids) {
        #warn "[$rid]\n";
        $s .= $bld->uid2row($rid)->rdb_row('data', $self->{'pad'},
                                           $self->{'gap'}) . "\n";
    }
    \$s;
}

###########################################################################
my $MSF_SEQ_WIDTH = 50;
my $MSF_PAD = '.';
my $MSF_GAP = '.';

#return alignment in MSF format
sub msf {
    my $self = shift;
    my ($bld, $aln, $s) = ($self->{'build'}, $self->{'align'}, '');

    my $CHECKSUM = '--------------------------------------&---*---.-----------------@ABCDEFGHIJKLMNOPQRSTUVWXYZ------ABCDEFGHIJKLMNOPQRSTUVWXYZ---~---------------------------------------------------------------------------------------------------------------------------------';

    my $checksum = sub {
	my $s = shift;
	my ($sum, $ch) = (0, 0);
	my $len = length($$s);
	while ($len--) {
	    $ch = ord substr($$s, $len, 1);
	    $ch = substr($CHECKSUM, $ch, 1);
	    $sum += (($len % 57) + 1) * ord $ch  if $ch ne '-';
	}
	$sum % 10000;
    };

    my $now = `date '+%B %d, %Y %H:%M'`;
    $now =~ s/\s0(\d{1})/ $1/; chomp $now; #padding %-d may not work

    if ($self->{'moltype'} eq 'aa') {
        $s .= "!!AA_MULTIPLE_ALIGNMENT 1.0\n";
    } else {
        $s .= "!!NA_MULTIPLE_ALIGNMENT 1.0\n";
    }
    $s .= "PileUp (MView)\n\n";
    $s .= sprintf("   MSF: %5d  Type: %s  %s  Check: %4d  ..\n\n",
          $aln->length, ($self->{'moltype'} eq 'aa' ? 'P' : 'N'), $now, 0);

    my $w = 0;
    foreach my $rid ($aln->visible_ids) {
        $w = length($bld->uid2row($rid)->cid)
            if length($bld->uid2row($rid)->cid) > $w;
    }

    foreach my $rid ($aln->visible_ids) {
        $s .=
            sprintf(" Name: %-${w}s Len: %5d  Check: %4d  Weight:  %4.2f\n",
                    $bld->uid2row($rid)->cid, $aln->length,
                    &$checksum(\$bld->uid2row($rid)->seq), 1.0);
    }
    $s .= "\n//\n\n";

    my ($pad, $gap, %seq) = ($self->{'pad'}, $self->{'gap'});

    foreach my $rid ($aln->visible_ids) {
        $seq{$rid} = $bld->uid2row($rid)->seq($pad, $gap);
    }

  LOOP:
    for (my $from = 0; ;$from += $MSF_SEQ_WIDTH) {
        my $ruler = 1;
        foreach my $rid ($aln->visible_ids) {
            last LOOP    if $from >= length($seq{$rid});
            my $tmp = substr($seq{$rid}, $from, $MSF_SEQ_WIDTH);
            my $tmplen = length($tmp);
            if ($ruler) {
                my $lo = $from + 1; my $hi = $from + $tmplen;
                $ruler = $tmplen - length("$lo") - length("$hi");
                $ruler = 1  if $ruler < 1;
                my $insert = int($tmplen / 10);
                $insert -= 1  if $tmplen % 10 == 0;
                $insert += $ruler;
                $insert = sprintf("%d%s%d", $lo, ' ' x $insert, $hi);
                $s .= sprintf("%-${w}s $insert\n", '');
                $ruler = 0;
            }
            $s .= sprintf("%-${w}s ", $bld->uid2row($rid)->cid);
            for (my $lo=0; $lo<$tmplen; $lo+=10) {
                $s .= substr($tmp, $lo, 10);
                $s .= ' '    if $lo < 40;
            }
            $s .= "\n";
        }
        $s .= "\n";
    }
    \$s;
}

###########################################################################
my $CLUSTAL_SEQ_WIDTH    = 60;
my $CLUSTAL_NAME_WIDTH   = 16;
my $CLUSTAL_RULER        = 1;
my $CLUSTAL_CONSERVATION = 1;

#return alignment in CLUSTAL/aln format
sub clustal {
    my $self = shift;
    my ($bld, $aln, $s) = ($self->{'build'}, $self->{'align'}, '');

    my $symcount = sub {
	my ($s, $pad, $gap, $c) = (@_, 0);
	my @s = split('', $$s);
	for (my $i=0; $i<@s; $i++) {
	    $c += 1  unless $s[$i] eq $pad or $s[$i] eq $gap;
	}
	$c;
    };

    $s .= "CLUSTAL 2.1 multiple sequence alignment (MView)\n\n\n";

    my $w = $CLUSTAL_NAME_WIDTH;
    my ($pad, $gap) = ($self->{'pad'}, $self->{'gap'});
    my @ids = $aln->visible_ids;
    my %seq = ();
    my %len = ();

    foreach my $rid (@ids) {
        #warn $rid, $bld->uid2row($rid), "\n";
        $w = length($bld->uid2row($rid)->cid)+1 if
             length($bld->uid2row($rid)->cid) >= $w;
        $seq{$rid} = $bld->uid2row($rid)->seq($pad, $gap);
        $len{$rid} = 0;
    }

  LOOP:
    for (my $from = 0; ;$from+=$CLUSTAL_SEQ_WIDTH) {
        foreach my $rid ($aln->visible_ids) {
            last LOOP  if $from >= length($seq{$rid});
            my $tmp = substr($seq{$rid}, $from, $CLUSTAL_SEQ_WIDTH);
            $s .= sprintf("%-${w}s", $bld->uid2row($rid)->cid);
            $s .= $tmp;
            if ($CLUSTAL_RULER) {
                my $syms = &$symcount(\$tmp, $pad, $gap);
                my $hi = $len{$rid} + $syms;
                $s .= sprintf(" %-d", $hi)  if $syms > 0;
                $len{$rid} = $hi;
            }
            $s .= "\n";
        }

        if ($CLUSTAL_CONSERVATION) {
            $s .= sprintf("%-${w}s", '');
            $s .= ${$aln->conservation(\@ids, $from + 1,
                                       $from + $CLUSTAL_SEQ_WIDTH,
                                       $self->{'moltype'})};
        }
        $s .= "\n\n";
    }
    \$s;
}


###########################################################################
1;
