# Makefile for Travis CI tests and installation
# https://docs.travis-ci.com/user/languages/perl

test:
	perl test_get_homologues.t
	rm -rf sample_plasmids_gbk_homologues/tmp/ 
	rm -rf sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algOMCL_e0_
	rm -f sample_plasmids_gbk_homologues/UnculturedbacteriumplasmidpRSB203_f0_0taxa_algOMCL_e0_.cluster_list
	rm -rf tmpclusters

install:
	perl install.pl

install_auto:
	perl install.pl no_databases
