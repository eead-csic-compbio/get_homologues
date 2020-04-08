# Copyright (C) 1997 Christophe Leroy, 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::SRS;

use Bio::MView::Option::Parameters;  #for $PAR

###########################################################################
#NCBI entrez 2007
my $URL_NCBI = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi';

my %MAP_NCBI = (
    #protein databases
    'dbj'     => "$URL_NCBI?val=&AC",
    'emb'     => "$URL_NCBI?val=&AC",
    'gb'      => "$URL_NCBI?val=&AC",
    'pdb'     => "$URL_NCBI?val=&AC&ID",
    'pir'     => "$URL_NCBI?val=&ID",
    'ref'     => "$URL_NCBI?val=&AC",
    'sp'      => "$URL_NCBI?val=&AC",
    );

###########################################################################
#EBI SRS server 2008
my $URL_EBI = 'http://srs.ebi.ac.uk/srsbin/cgi-bin/wgetz';

my %MAP_EBI = (
    #protein databases
    'uniprot'   => "$URL_EBI?-e+[uniprot:(&ID)]|[uniprot:(&AC)]+-noSession",
    'swiss'     => "$URL_EBI?-e+[uniprot:(&ID)]|[uniprot:(&AC)]+-noSession",
    'sw'        => "$URL_EBI?-e+[uniprot:(&ID)]|[uniprot:(&AC)]+-noSession",
    'sp'        => "$URL_EBI?-e+[uniprot:(&ID)]|[uniprot:(&AC)]+-noSession",
    'uniref100' => "$URL_EBI?-e+[uniref100:(&AC)]+-noSession",
    'uniref90'  => "$URL_EBI?-e+[uniref90:(&AC)]+-noSession",
    'uniref50'  => "$URL_EBI?-e+[uniref50:(&AC)]+-noSession",
    'ur100'     => "$URL_EBI?-e+[uniref100:(&AC)]+-noSession",
    'ur90'      => "$URL_EBI?-e+[uniref90:(&AC)]+-noSession",
    'ur50'      => "$URL_EBI?-e+[uniref50:(&AC)]+-noSession",
    'uniparc'   => "$URL_EBI?-e+[uniparc:(&ID)]|[uniparc:(&AC)]+-noSession",
    'ipi'       => "$URL_EBI?-e+([ipi:(&ID)]|[ipi:(&AC)])+-noSession",
    'pdb'       => "$URL_EBI?-e+[pdb:(&ID)]+-noSession",

    #nucleotide databases
    'EM_REL' => "$URL_EBI?-e+([emblidacc:(&ID)]%3Eembl)|[embl:(&AC)]+-view+EmblEntry+-noSession",
    'EM_NEW' => "$URL_EBI?-e+([emblidacc:(&ID)]%3Eembl)|[embl:(&AC)]+-view+EmblEntry+-noSession",
    );

#EMBL nucleotide database aliases to EM_REL
$MAP_EBI{'EM_ANN'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_CON'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_EST'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_GSS'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_HTC'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_HTG'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_PAT'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_STD'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_STS'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_TPA'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_WGS'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_ENV'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_FUN'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_HUM'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_INV'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_MAM'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_MUS'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_PHG'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_PLN'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_PRO'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_ROD'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_SYN'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_TGN'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_UNC'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_VRL'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_VRT'} = \$MAP_EBI{'EM_REL'};
$MAP_EBI{'EM_OM'}  = \$MAP_EBI{'EM_REL'};

######################################################################
# public methods
######################################################################
sub srsLink {
    my ($tag) = @_;
    return ''  unless $PAR->get('srs');

    my $link = '';

    #EDIT: PLACE SITE SPECIFIC PATTERNS FIRST

    #MSKCC style names: db|type|val
    if ($tag =~ /^[^|:]+\|(?:id|ac)\|/) {
        if ((my @tmp = split(/\|/, $tag)) >= 2) {
            my ($db, $type, $val) = ($tmp[0], $tmp[1], $tmp[2]);
            $link = getLink(\%MAP_EBI, $db, $val, $val);
        }
        #warn "Split MSKGG $tag => ($db, $type, $val) => $link\n";
        return $link;
    }

    #EDIT: SITE SPECIFIC PATTERNS END

    #NCBI style names: gi|number|db|ac|id
    if ($tag =~ /^gi\|[0-9]+\|/) {
        if ((my @tmp = split(/\|/, $tag)) >= 5) {
            my ($db, $ac, $id) = ($tmp[2], $tmp[3], $tmp[4]);
            $link = getLink(\%MAP_NCBI, $db, $ac, $id);
        }
        #warn "Split NCBI (long) $tag => ($db, $ac, $id) => $link\n";
        return $link;
    }

    #NCBI style names: db|ac|id,  db|ac|
    if ($tag =~ /^[^|:]+\|/) {
        if ((my @tmp = split(/\|/, $tag)) >= 2) {
            my ($db, $ac, $id) = ($tmp[0], $tmp[1], $tmp[2]||"");
            if ($ac =~ /^(\S+):\S+/) {
                #SEGS style name:range
                $ac = $id = $1;
            }
            $link = getLink(\%MAP_NCBI, $db, $ac, $id);
        }
        #warn "Split NCBI (short) $tag => ($db, $ac, $id) => $link\n";
        return $link;
    }

    #EBI/GCG style names: db:id
    if ($tag =~ /^[^|:]+\:/) {
        if ((my @tmp = split(/:/, $tag)) >= 2) {
            my ($db, $ac, $id) = ($tmp[0], $tmp[1], $tmp[1]);
            $link = getLink(\%MAP_EBI, $db, $ac, $id);
        }
        #warn "Split EBI $tag => ($db,$ac,$id) => $link\n";
        return $link;
    }

    warn "No SRS match for '$tag'";

    #no match
    return $link;
}

######################################################################
# private methods
######################################################################
sub getLink() {
    my ($map, $db, $ac, $id) = @_;
    my $link = '';

    foreach $db ($db, lc $db, uc $db) {
        if (exists $map->{$db}) {
            $link = $map->{$db};
            last;
        }
    }

    if (ref $link) { #follow aliases
        if (defined $$link) {
            $link = $$link;
        } else {
            $link = '';
        }
    }

    #remove trailing punctuation
    if ($ac =~ /^(.+[a-zA-Z0-9])[-.,;:_|]+$/) {
        $ac = $1;
    }
    if ($id =~ /^(.+[a-zA-Z0-9])[-.,;:_|]+$/) {
        $id = $1;
    }

    #remove trailing single character qualifiers: .1 _A
    if ($ac =~ /^([^.]+)[._].$/) {
        $ac = $1;
    }
    if ($id =~ /^([^.]+)[._].$/) {
        $id = $1;
    }

    #substitute into template
    $link =~ s/\&AC/$ac/g;
    $link =~ s/\&ID/$id/g;
    #warn "Link  ($db,$ac,$id) => $link\n";

    return $link;
}

###########################################################################
1;
