# Makefile for tests and installation

test:
	perl test_get_homologues.t

testR:
	perl test_get_homologues.t testR

test_swiss:
	perl test_get_homologues.t swiss

test_nonet:
	perl test_get_homologues.t nonet

clean:	
	rm -rf sample_plasmids_gbk_homologues/tmp/ 
	rm -rf sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algOMCL_e0_
	rm -f sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algOMCL_e0_.cluster_list
	rm -rf sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algCOG_e0_
	rm -f sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algCOG_e0_.cluster_list
	rm -f sample_intersection/*.png
	rm -f _sample_cluster.fna_* sample_cluster.fna_*
	rm -rf sample_transcripts_fasta_est_homologues/

cleanR:
	rm -f pangenome_matrix_variable_sites_only.tsv gower_dist_matrix.tab hclust_gower-ward* hcluster_ward*
	rm -f height_* hcluster_pangenome_matrix_script_run_at_* ANDg_* plot_matrix_heatmap_script_run_at*

install:
	perl install.pl

install_swissprot:
	perl install.pl swissprot

install_auto:
	perl install.pl no_databases
