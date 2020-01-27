use strict;
use warnings;
use Test::More tests => 28;
use FindBin '$Bin';
use lib "$Bin/../lib";
use lib "$Bin/../lib/bioperl-1.5.2_102/";

BEGIN { use_ok('ForkManager') };

BEGIN { use_ok('HPCluster') };

BEGIN { use_ok('phyTools') };

BEGIN { use_ok('marfil_homology') };

BEGIN { use_ok('marfil_homology') };

BEGIN { use_ok('marfil_homology') };

ok( eval{ `perl $Bin/../get_homologues.pl -v ` } =~ /Checking required binaries/ , 'get_homologues.pl -v' );

ok( eval{ `perl $Bin/../get_homologues-est.pl -v ` } =~ /Checking required binaries/ , 'get_homologues-est.pl -v' );

ok( eval{ `perl $Bin/../add_pancore_matrices.pl ` } =~ /\[options\]/ , 'add_pancore_matrices.pl' );

ok( eval{ `perl $Bin/../add_pangenome_matrices.pl` } =~ /\[options\]/ , 'add_pangenome_matrices.pl' );

ok( eval{ `perl $Bin/../check_BDBHs.pl` } =~ /this message/ , 'check_BDBHs.pl' );

ok( eval{ `perl $Bin/../_cluster_makeHomolog.pl` } =~ /\[options\]/ , '_cluster_makeHomolog.pl' );

ok( eval{ `perl $Bin/../_cluster_makeInparalog.pl` } =~ /\[options\]/ , '_cluster_makeInparalog.pl' );

ok( eval{ `perl $Bin/../_cluster_makeIsoform.pl` } =~ /\[options\]/ , '_cluster_makeIsoform.pl' );

ok( eval{ `perl $Bin/../_cluster_makeOrtholog.pl` } =~ /\[options\]/ , '_cluster_makeOrtholog.pl' );

ok( eval{ `perl $Bin/../compare_clusters.pl` } =~ /\[options\]/ , 'compare_clusters.pl' );

ok( eval{ `perl $Bin/../download_genomes_ncbi.pl test` } =~ /\$NCBIHOST ok/ , 'download_genomes_ncbi.pl test' );

ok( eval{ `$Bin/../hcluster_pangenome_matrix.sh` } =~ /check_dependencies/ , 'hcluster_pangenome_matrix.sh' );

ok( eval{ `perl $Bin/../install.pl test` } =~ /testing only/ , 'install.pl test' );

ok( eval{ `perl $Bin/../make_nr_pangenome_matrix.pl` } =~ /\[options\]/ , 'make_nr_pangenome_matrix.pl' );

ok( eval{ `perl $Bin/../parse_pangenome_matrix.pl` } =~ /\[options\]/ , 'parse_pangenome_matrix.pl' );

ok( eval{ `perl $Bin/../pfam_enrich.pl` } =~ /\[options\]/ , 'pfam_enrich.pl' );

ok( eval{ `$Bin/../plot_matrix_heatmap.sh` } =~ /check_dependencies/ , 'plot_matrix_heatmap.sh' );

ok( eval{ `perl $Bin/../plot_pancore_matrix.pl` } =~ /\[options\]/ , 'plot_pancore_matrix.pl' );

ok( eval{ `perl $Bin/../_split_blast.pl` } =~ /Usage/ , '_split_blast.pl' );

ok( eval{ `perl $Bin/../_split_hmmscan.pl` } =~ /Usage/ , '_split_hmmscan.pl' );

ok( eval{ `perl $Bin/../transcripts2cds.pl` } =~ /\[options\]/ , 'transcripts2cds.pl' );
# this implicitly  tests modules in libs/est/

ok( eval{ `perl $Bin/../transcripts2cdsCPP.pl` } =~ /\[options\]/ , 'transcripts2cdsCPP.pl' );
# this needs optional module Inline::CPP
