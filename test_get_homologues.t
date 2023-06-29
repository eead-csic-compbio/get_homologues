#!/usr/bin/env perl 

use strict;
use warnings;
use Test::More;
use FindBin '$Bin';
use lib "$Bin/lib";
use lib "$Bin/lib/bioperl-1.5.2_102/";

chdir($Bin);

BEGIN { use_ok('ForkManager') };

BEGIN { use_ok('HPCluster') };

BEGIN { use_ok('phyTools') };

BEGIN { use_ok('marfil_homology') };

my ($testR,$testSP) = (0,0);
if(defined($ARGV[0])){
	if($ARGV[0] eq 'testR') { 
		$testR = 1 
	}
	elsif($ARGV[0] eq 'swiss'){ 
		$testSP = 1;
	}
}

ok( eval{ `perl ./add_pancore_matrices.pl ` } =~ /\[options\]/ , 'add_pancore_matrices.pl' );

ok( eval{ `perl ./add_pangenome_matrices.pl` } =~ /\[options\]/ , 'add_pangenome_matrices.pl' );

ok( eval{ `perl ./check_BDBHs.pl` } =~ /this message/ , 'check_BDBHs.pl' );

ok( eval{ `perl ./_cluster_makeHomolog.pl` } =~ /\[options\]/ , '_cluster_makeHomolog.pl' );

ok( eval{ `perl ./_cluster_makeInparalog.pl` } =~ /\[options\]/ , '_cluster_makeInparalog.pl' );

ok( eval{ `perl ./_cluster_makeIsoform.pl` } =~ /\[options\]/ , '_cluster_makeIsoform.pl' );

ok( eval{ `perl ./_cluster_makeOrtholog.pl` } =~ /\[options\]/ , '_cluster_makeOrtholog.pl' );

if($ARGV[0] ne 'nonet') {
  ok( eval{ `perl ./download_genomes_ncbi.pl test` } =~ /\$NCBIHOST ok/ , 'download_genomes_ncbi.pl test' );
} else {
  ok( eval{ `perl ./download_genomes_ncbi.pl 2>&1 ` } =~ /usage/ , 'download_genomes_ncbi.pl' );  
}

ok( eval{ `perl ./install.pl test` } =~ /testing only/ , 'install.pl test' );

ok( eval{ `perl ./make_nr_pangenome_matrix.pl` } =~ /\[options\]/ , 'make_nr_pangenome_matrix.pl' );

ok( eval{ `perl ./pfam_enrich.pl` } =~ /\[options\]/ , 'pfam_enrich.pl' );

ok( eval{ `perl ./plot_pancore_matrix.pl` } =~ /\[options\]/ , 'plot_pancore_matrix.pl' );

ok( eval{ `perl ./_split_blast.pl` } =~ /Usage/ , '_split_blast.pl' );

ok( eval{ `perl ./_split_hmmscan.pl` } =~ /Usage/ , '_split_hmmscan.pl' );

ok( eval{ `perl ./transcripts2cds.pl` } =~ /\[options\]/ , 'transcripts2cds.pl' );
# implicitly tests modules in libs/est/

if($testSP){
	# requires optional module Inline::CPP, not available in bioconda May2022
	ok( eval{ `perl ./transcripts2cdsCPP.pl sample_cluster.fna` } =~ /sample_cluster.fna_l50_E1e-05.cds/ , 'transcripts2cdsCPP.pl' );
} 

ok( eval{ `perl ./get_homologues-est.pl -v ` } =~ /Checking required binaries/ , 'get_homologues-est.pl -v' );

ok( eval{ `perl ./get_homologues-est.pl -d sample_transcripts_fasta -m dryrun` } =~ /check the list of pending commands/ , 'get_homologues-est.pl -dryrun' );

ok( eval{ `perl ./get_homologues.pl -d sample_plasmids_gbk -M -t 0` } =~ /number_of_clusters = 19\d+/, 'get_homologues.pl -d sample_plasmids_gbk -M' );

ok( eval{ `perl ./get_homologues.pl -d sample_plasmids_gbk -G -t 0` } =~ /number_of_clusters = 19\d+/, 'get_homologues.pl -d sample_plasmids_gbk -G' );

ok( eval{ `perl ./compare_clusters.pl -d sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algOMCL_e0_,sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algCOG_e0_ -o sample_intersection -m` } =~ /pangenome_matrix_t0.tab/ , 'compare_clusters.pl ' );

ok( eval{ `perl ./parse_pangenome_matrix.pl -m sample_intersection/pangenome_matrix_t0.tab -s` } =~ /sample_intersection\/pangenome_matrix_t0__shell.png/, 'parse_pangenome_matrix.pl -s' );

# requires module GD and package libgd-dev
ok( eval{ `perl ./parse_pangenome_matrix.pl -m sample_intersection/pangenome_matrix_t0.tab -A sample_plasmids_gbk/A.txt -B sample_plasmids_gbk/B.txt -g -p 'Klebsiella oxytoca KOX105'` } =~ /chromosome\/contig/, 'parse_pangenome_matrix.pl -p' );

# uses lib/mview
ok( eval{ `perl ./annotate_cluster.pl -f sample_cluster.fna -b 2>&1` } =~ /taxa included in alignment: 14/, 'annotate_cluster.pl -f sample_cluster -b' );

# requires R dependencies
if($testR == 1){
	ok( eval{ `bash ./hcluster_pangenome_matrix.sh -i sample_intersection/pangenome_matrix_t0.tab 2>&1` } =~ 
		/heatmap.pdf was generated/ , 'hcluster_pangenome_matrix.sh -i ' );
} else {
	ok( eval{ `bash ./hcluster_pangenome_matrix.sh` } =~ /synopsis/ , 'hcluster_pangenome_matrix.sh' );
}

if($testR == 1){
	ok( eval{ `bash ./plot_matrix_heatmap.sh -i sample_intersection/pangenome_matrix_t0.tab -o pdf 2>&1` } =~ 
	/heatmap.pdf was produced/ , 'plot_matrix_heatmap.sh -i' );
} else {
	ok( eval{ `bash ./plot_matrix_heatmap.sh` } =~ /synopsis/ , 'plot_matrix_heatmap.sh' );
}

if($testSP){
	done_testing(30)
} else {
	done_testing(29)
}
