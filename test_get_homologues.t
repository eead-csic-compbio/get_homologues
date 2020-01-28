use strict;
use warnings;
use Test::More tests => 28;
use lib "lib";
use lib "lib/bioperl-1.5.2_102/";

BEGIN { use_ok('ForkManager') };

BEGIN { use_ok('HPCluster') };

BEGIN { use_ok('phyTools') };

BEGIN { use_ok('marfil_homology') };

BEGIN { use_ok('marfil_homology') };

BEGIN { use_ok('marfil_homology') };

ok( eval{ `perl ./add_pancore_matrices.pl ` } =~ /\[options\]/ , 'add_pancore_matrices.pl' );

ok( eval{ `perl ./add_pangenome_matrices.pl` } =~ /\[options\]/ , 'add_pangenome_matrices.pl' );

ok( eval{ `perl ./check_BDBHs.pl` } =~ /this message/ , 'check_BDBHs.pl' );

ok( eval{ `perl ./_cluster_makeHomolog.pl` } =~ /\[options\]/ , '_cluster_makeHomolog.pl' );

ok( eval{ `perl ./_cluster_makeInparalog.pl` } =~ /\[options\]/ , '_cluster_makeInparalog.pl' );

ok( eval{ `perl ./_cluster_makeIsoform.pl` } =~ /\[options\]/ , '_cluster_makeIsoform.pl' );

ok( eval{ `perl ./_cluster_makeOrtholog.pl` } =~ /\[options\]/ , '_cluster_makeOrtholog.pl' );

ok( eval{ `perl ./download_genomes_ncbi.pl test` } =~ /\$NCBIHOST ok/ , 'download_genomes_ncbi.pl test' );

ok( eval{ `bash ./hcluster_pangenome_matrix.sh` } =~ /check_dependencies/ , 'hcluster_pangenome_matrix.sh' );

ok( eval{ `perl ./install.pl test` } =~ /testing only/ , 'install.pl test' );

ok( eval{ `perl ./make_nr_pangenome_matrix.pl` } =~ /\[options\]/ , 'make_nr_pangenome_matrix.pl' );

ok( eval{ `perl ./pfam_enrich.pl` } =~ /\[options\]/ , 'pfam_enrich.pl' );

ok( eval{ `bash ./plot_matrix_heatmap.sh` } =~ /check_dependencies/ , 'plot_matrix_heatmap.sh' );

ok( eval{ `perl ./plot_pancore_matrix.pl` } =~ /\[options\]/ , 'plot_pancore_matrix.pl' );

ok( eval{ `perl ./_split_blast.pl` } =~ /Usage/ , '_split_blast.pl' );

ok( eval{ `perl ./_split_hmmscan.pl` } =~ /Usage/ , '_split_hmmscan.pl' );

ok( eval{ `perl ./transcripts2cds.pl` } =~ /\[options\]/ , 'transcripts2cds.pl' );
# implicitly tests modules in libs/est/

ok( eval{ `perl ./transcripts2cdsCPP.pl` } =~ /\[options\]/ , 'transcripts2cdsCPP.pl' );
# requires optional module Inline::CPP

ok( eval{ `perl ./get_homologues-est.pl -v ` } =~ /Checking required binaries/ , 'get_homologues-est.pl -v' );

ok( eval{ `perl ./get_homologues.pl -d sample_plasmids_gbk -M -t 0` } =~ /number_of_clusters = 19\d+/, 'get_homologues.pl -d sample_plasmids_gbk -M' );

ok( eval{ `perl ./compare_clusters.pl -d sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algOMCL_e0_ -o tmpclusters -m` } =~ /pangenome_matrix_t0.tab/ , 'compare_clusters.pl -d sample_plasmids_gbk_homologues/...' );

ok( eval{ `perl ./parse_pangenome_matrix.pl -m tmpclusters/pangenome_matrix_t0.tab -s` } =~ /pan-genome size estimates/, 'parse_pangenome_matrix.pl' );
