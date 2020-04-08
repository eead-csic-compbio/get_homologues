# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
#
# NCBI BLAST 2.0, PSI-BLAST 2.0, BLAST+
#
#   blastp, blastn, blastx, tblastn, tblastx
#
###########################################################################
package Bio::MView::Build::Format::BLAST2;

use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Build::Format::BLAST0;

use strict;
use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST0);

sub begin_parse {
    my $self = shift;
    #warn "blast2::begin_parse: cycle ", $self->cycle,
    #    " strand ", $self->strand, " schedule ", $self->{scheduler}->item, "\n";
    my $s = $self->{'entry'}->parse("SEARCH", $self->cycle);
    return ()  unless defined $s;
    my $r = $s->parse(qw(RANK));
    return ()  unless defined $r;
    my $h  = $self->{'entry'}->parse(qw(HEADER));
    ($s, $h, $r);
}

sub end_parse {
    my $self = shift;
    #warn "blast2::end_parse: free SEARCH\n";
    $self->{'entry'}->free_parsers(qw(SEARCH RANK HEADER));
}

sub parse_record {
    my ($self, $k, $hdr, $sum, $aln) = @_;

    #sub str { join ", ", map { defined $_ ? $_ : 'undef' } @_ }
    #warn "[@{[str($k, $hdr, $sum, $aln)]}]\n";

    sub load_hdr_extra {
        my ($hdr, $info) = @_;
        my $keep = {};
        for (my $i=0; $i < @{$hdr->{'fields'}}; $i++) {
            my $k = $hdr->{'fields'}->[$i];
            next  unless exists $hdr->{'extra'}->{$k};  #unwanted field
            next  unless $hdr->{'counts'}->[$i] > 0;    #zero occupancy
            $keep->{$k} = 1;
            $info->{$k} = $hdr->{'extra'}->{$k};
        }
        return $keep;
    }

    sub load_row_extra {
        my ($keep, $row, $info) = @_;
        while (my ($k, $v) = each %{$row->{'extra'}}) {
            next  unless exists $keep->{$k};
            $info->{$k} = $v;
        }
    }

    my ($id, $desc, $bits, $expect, $n) = ('','','','','');

    $k = ''  unless defined $k;  #row number = rank

    #not mutually exclusive: order of preference
    if (defined $hdr) {
        ($id, $desc) = ($hdr->{'query'}, $hdr->{'summary'});
    } elsif (defined $sum) {
        ($id, $desc) = ($sum->{'id'}, $sum->{'desc'});
    } elsif (defined $aln) {
        ($id, $desc) = ($aln->{'id'}, $aln->{'summary'});
    }

    #always
    if (defined $aln) {
        ($bits, $expect, $n) = ($aln->{'bits'}, $aln->{'expect'}, $aln->{'n'});
    }

    my $info1 = {
        'cycle'  => $self->cycle,
        'bits'   => $bits,
        'expect' => $expect,
        'n'      => $n,
    };

    my $info2 = {};

    #extract any optional fields from BLAST_OF7 skipping any fields
    #with zero occupancy
    if (defined $hdr and exists $hdr->{'extra'}) {
        $self->{'extra_keys'} = load_hdr_extra($hdr, $info2);
    } elsif (defined $aln and exists $aln->{'extra'}) {
        load_row_extra($self->{'extra_keys'}, $aln, $info2);
    }
    #use Bio::Util::Object qw(dump_hash); warn dump_hash($info);

    return ($k, $id, $desc, $info1, $info2);
}

sub get_scores {
    my ($self, $list) = @_;
    my @tmp = $self->SUPER::get_scores($list);
    return ()  unless @tmp;

    my $bits = \$tmp[1];

    #try to preserve original precision/format as far as possible
    my $val = sprintf("%.0f", $$bits);     #bits as integer
    if ($val != $$bits) {                  #num. different? put 1dp back
        $$bits = sprintf("%.1f", $$bits);  #bit score to 1 dp
    }

    @tmp;
}

#score/significance filter
sub skip_hsp {
    my ($self, $hsp) = @_;
    return 1  if
        defined $PAR->get('minbits') and $hsp->{'bits'} < $PAR->get('minbits');
    return 1  if
        defined $PAR->get('maxeval') and $hsp->{'expect'} > $PAR->get('maxeval');
    return 0;
}

#standard attribute names
sub attr_score { 'bits' }
sub attr_sig   { 'expect' }


###########################################################################
package Bio::MView::Build::Row::BLAST2;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::BLAST);

#suppress cycle in all output; it's only for blastp/psi-blastp
sub schema {[
    # use? rdb?  key              label          format   default
    [ 0,   0,    'cycle',         'cycle',       '2N',    ''   ],
    [ 2,   2,    'bits',          'bits',        '5N',    ''   ],
    [ 3,   3,    'expect',        'E-value',     '9N',    ''   ],
    [ 4,   4,    'n',             'N',           '2N',    ''   ],
    [ 5,   5,    'query_orient',  'qy',          '2S',    '?'  ],
    [ 6,   6,    'sbjct_orient',  'ht',          '2S',    '?'  ],
    ]
}

sub optional_schema {[
    # use? rdb?  key              label          format   default
    [ 10,  10,   'staxid',        'staxid',      '6S',    ''   ],
    [ 11,  11,   'ssciname',      'ssciname',    '8S',    ''   ],
    [ 12,  12,   'scomname',      'scomname',    '8S',    ''   ],
    [ 13,  13,   'sblastname',    'sblastname',  '10S',   ''   ],
    [ 14,  14,   'sskingdom',     'sskingdom',   '9S',    ''   ],
    [ 15,  15,   'staxids',       'staxids',     '7S',    ''   ],
    [ 16,  16,   'sscinames',     'sscinames',   '9S',    ''   ],
    [ 17,  17,   'scomnames',     'scomnames',   '9S',    ''   ],
    [ 18,  18,   'sblastnames',   'sblastnames', '11S',   ''   ],
    [ 19,  19,   'sskingdoms',    'sskingdoms',  '10S',   ''   ],
    ]
}


###########################################################################
package Bio::MView::Build::Row::BLAST2::blastp;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::BLAST2);

#enable cycle in rdb tabulated output
#suppress query and sbjct orientations for blastp data
sub schema {[
    # use? rdb?  key              label          format   default
    [ 0,   1,    'cycle',         'cycle',       '2N',    ''   ],
    [ 2,   2,    'bits',          'bits',        '5N',    ''   ],
    [ 3,   3,    'expect',        'E-value',     '9N',    ''   ],
    [ 4,   4,    'n',             'N',           '2N',    ''   ],
    [ 0,   0,    'query_orient',  'qy',          '2S',    '?'  ],
    [ 0,   0,    'sbjct_orient',  'ht',          '2S',    '?'  ],
    ]
}


###########################################################################
package Bio::MView::Build::Row::BLAST2::blastn;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::BLAST2);


###########################################################################
package Bio::MView::Build::Row::BLAST2::blastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::BLASTX);

sub schema { Bio::MView::Build::Row::BLAST2::schema }
sub optional_schema { Bio::MView::Build::Row::BLAST2::optional_schema }


###########################################################################
package Bio::MView::Build::Row::BLAST2::tblastn;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::BLAST2);


###########################################################################
package Bio::MView::Build::Row::BLAST2::tblastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Row::BLASTX);

sub schema { Bio::MView::Build::Row::BLAST2::schema }
sub optional_schema { Bio::MView::Build::Row::BLAST2::optional_schema }


###########################################################################
###########################################################################
package Bio::MView::Build::Format::BLAST2::blastp;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST2
          Bio::MView::Build::Format::BLAST0::blastp
);

sub scheduler { 'cycle' }

sub record_type { 'Bio::MView::Build::Row::BLAST2::blastp' }

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s  if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Search cycle: " . $self->cycle . "\n";
    $s;
}


###########################################################################
package Bio::MView::Build::Format::BLAST2::blastn;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST2
          Bio::MView::Build::Format::BLAST0::blastn
);

sub scheduler { 'strand' }

sub record_type { 'Bio::MView::Build::Row::BLAST2::blastn' }

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s  if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Query orientation: " . $self->strand . "\n";
    $s;
}


###########################################################################
package Bio::MView::Build::Format::BLAST2::blastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST2
          Bio::MView::Build::Format::BLAST0::blastx
);

sub scheduler { 'strand' }

sub record_type { 'Bio::MView::Build::Row::BLAST2::blastx' }

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s    if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Query orientation: " . $self->strand . "\n";
    $s;
}


###########################################################################
package Bio::MView::Build::Format::BLAST2::tblastn;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST2
          Bio::MView::Build::Format::BLAST0::tblastn
);

sub scheduler { 'none' }

sub record_type { 'Bio::MView::Build::Row::BLAST2::tblastn' }


###########################################################################
package Bio::MView::Build::Format::BLAST2::tblastx;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Build::Format::BLAST2
          Bio::MView::Build::Format::BLAST0::tblastx
);

sub scheduler { 'strand' }

sub record_type { 'Bio::MView::Build::Row::BLAST2::tblastx' }

sub subheader {
    my ($self, $quiet) = (@_, 0);
    my $s = '';
    return $s  if $quiet;
    $s  = $self->SUPER::subheader($quiet);
    $s .= "Query orientation: " . $self->strand . "\n";
    $s;
}


###########################################################################
1;
