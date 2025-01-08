package marfil_homology;

# version 2.1

# Code library created by Bruno Contreras-Moreira and Pablo Vinuesa (2006-2025)
# mainly for get_homologues.pl and get_homologues-est.pl

# Contains code originally part of the following software:

## ORTHOMCL [2007-04-04] Version 1.4 Li Li, Feng Chen <fengchen@sas.upenn.edu>
## Copyright (C) 2004~2006 by University of Pennsylvania, Philadelphia, PA USA.

## COGtriangles version 2.1
## Kristensen DM, Kannan L, Coleman MK, Wolf YI, Sorokin A, Koonin EV, Mushegian A. Bioinformatics 2010.

use strict;
require Exporter;
use phyTools;
eval{ require Storable } || die "marfil_homology: failed to import perl module Storable, please install it (please check manual)\n@_";
use Storable qw(retrieve nstore);
use File::Temp qw(tempfile);
use Symbol qw(gensym);
use sort 'stable';
use List::Util 'first';

our @ISA = qw( Exporter );

our @EXPORT = qw(
  executeFORMATDB format_BLASTP_command executeMCL constructDirectory constructAllFasta matrix
  blast_parse blast_parse_COG blastqueryab open_bpofile numeric_value simspan add_unmatched_singletons
  same_length makeInparalog matchlen makeParalog makeOrtholog makeHomolog findAllOrthologiesORTHMCL
  find_PARANOID_clusters simspan_hsps getline_from_bpofile construct_indexes construct_taxa_indexes
  find_local_alignment_coords save_params check_different_params parse_MCL_clusters
  pfam_parse construct_Pfam_hash cluster_lineage_expansions find_COGs sort_blast_results
  find_OMCL_clusters split_Pfam_clusters sort_taxon_hits format_BLASTN_command rename_blast_outfiles
  get_makeOrtholog_outfilename get_makeInparalog_outfilename get_makeHomolog_outfilename
  format_SPLITBLAST_command format_SPLITHMMPFAM_command format_HMMPFAM_command
  get_makeIsoform_outfilename makeIsoform nr_blast_report flag_small_clusters 
  construct_redundant_hash write_redundant_hash parse_Pfam_freqs get_string_with_previous_genomes
  format_BLASTP_command_aligns format_BLASTN_command_aligns 
  executeDIAMONDMAKEDB format_DIAMONDblastp_command 

  $BLAST_PVALUE_CUTOFF_DEFAULT
  $PERCENT_IDENTITY_CUTOFF_DEFAULT
  $PERCENT_IDENTITY_CUTOFF_EST_DEFAULT
  $PERCENT_MATCH_CUTOFF_DEFAULT
  $MIN_PERCENT_LENGTH_DIFFERENCE
  $MIN_PERSEQID_HOM
  $MIN_COVERAGE_HOM
  $MIN_PERSEQID_HOM_EST
  $MIN_COVERAGE_HOM_EST
  $SOFTCOREFRACTION
  $FASTAEXTENSION
  $BLAST_NOCPU
  $BLASTP
  $BLASTN
  $SPLITBLAST
  $HMMPFAM
  $SORTBIN

  $MAX_WEIGHT_DEFAULT
  $MCL_INFLATION_DEFAULT
  $MIN_BITSCORE_SIM_DEFAULT
  $MIN_NEIGHBOR_CORR_DEFAULT
  $MIN_BITSCORE_DEFAULT
  $MIN_SEGMENT_COVER_DEFAULT
  $lockfile
  $all_fa_file
  $selected_genomes_file
  $blast_file
  $bpo_file
  $blast_file_nr
  $bpo_file_nr
  $taxa_index_file
  $pfam_file
  $redundant_file
  $matrix_file
  $mcl_file
  $parameter_BDBH_log
  $parameter_INP_log
  $parameter_HOM_log
  $parameter_OMCL_log
  $parameter_PRND_log
  $parameter_COGS_log
  $parameter_ISOS_log
  $paranoid_file
  @taxa
  %gindex
  @gindex2
  %blastquery
  @graph
  %weight
  %taxa_bpo_index
  %pfam_hash
  %redundant
  $p2ofilename
  $lsefilename
  $lseoutfilename
  $coglogfilename
  $cogblasthits
);

# executable software
our $SORTLIMITRAM        = "500M"; # buffer size should fit in all computers
our $SORTBIN             = "sort --buffer-size=$SORTLIMITRAM"; # sort is part of all Linux systems, otherwise edit
our $GZIPBIN             = "gzip";

# all these defined in phyTools.pm:
our $BLASTP              = $ENV{"EXE_BLASTP"};
our $BLASTN              = $ENV{"EXE_BLASTN"};
our $FORMATDB            = $ENV{"EXE_FORMATDB"};
our $BLAST_NOCPU         = 2; # most current cpus are multicore
our $SPLITBLAST          = $ENV{"EXE_SPLITBLAST"};

our $DIAMONDMAKEDB       = $ENV{"EXE_DMNFT"};
our $DIAMONDPEXE         = $ENV{"EXE_DMNDP"};  

our $SPLITHMMPFAM        = $ENV{"EXE_SPLITHMMPFAM"};
our $HMMPFAM             = $ENV{"EXE_HMMPFAM"};
our $MCL                 = $ENV{"EXE_MCL"};
our $MULTIPARANOID       = $ENV{"EXE_MULTIPARANOID"};
our $MAKEHASHEXE         = $ENV{"EXE_MAKEHASH"};
our $READBLASTEXE        = $ENV{"EXE_READBLAST"};
our $LSEEXE              = $ENV{"EXE_COGLSE"};
our $COGTRIANGLESEXE     = $ENV{"EXE_COGTRI"};

# runtime & cutoff variables
our $BLAST_PVALUE_CUTOFF_DEFAULT         = 1e-5;
our $BLAST_DB_SIZE                       = 100000000; # as in COGs, to ensure comparable BLAST scores among diff runs
our $FASTAEXTENSION                      = 'fasta';
our $PERCENT_IDENTITY_CUTOFF_DEFAULT     = 1;         # Both PercentIdentity and PercentMatch cutoff are
our $PERCENT_MATCH_CUTOFF_DEFAULT        = 75.0;      # as in COGS, 70.0 for prok set to 0 by default in OrthoMCL, 0-100%, 50 used in INPARANOID
our $PERCENT_IDENTITY_CUTOFF_EST_DEFAULT = 95;        # for transcripts of same species, as in Trinity 
                                                      # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3571712

our $MIN_PERSEQID_HOM                    = 0.0;       # by default there's no identity cutoff, as long as BLAST returns a hit
our $MIN_COVERAGE_HOM                    = 20.0;      # by default a coverage of 20% of the longest sequence is enough to call homology
                                                      # Tettelin et al (2005): minimum %sequence identity required for homology = 50%
                                                      # Tettelin et al (2005): minimum %sequence coverage required for homology = 50%

our $MIN_PERSEQID_HOM_EST                = 70.0;      # The following figures correspond to BLASTN::megablast alignments: 
                                                      # 5th-percentil of barley transcripts = 84.62
                                                      # 5th-percentil of A.thaliana CDS = 78.77
our $MIN_COVERAGE_HOM_EST                = 50.0;      # This coverage allows matching sequences with retained/unprocessed introns

our $MIN_PERCENT_LENGTH_DIFFERENCE       = 80.0;      # 70.0 for prok short/long 0-100% [not benchmarked]; not used in INPARANOID 
our $MCL_INFLATION_DEFAULT               = 1.5;       # papers mention values from 1 to 4
our $MAX_WEIGHT_DEFAULT                  = 316;       # minimum e-value of your blast; ie if 9e-99 you should use -log(9e-99)=100.
our $MIN_BITSCORE_SIM_DEFAULT            = 28;        # Song,Joseph,Davis,Durand(2008)PlOS PMID:18475320
our $MIN_NEIGHBOR_CORR_DEFAULT           = 0.0;       # 0.4 in Song,Joseph,Davis,Durand(2008)PlOS PMID:18475320
our $MIN_BITSCORE_DEFAULT                = 50;        # in bits, as in inparanoid 3.0
our $MIN_SEGMENT_COVER_DEFAULT           = 75;        # default % in COGtriangles, (50 in inparanoid 3.0)
our $SOFTCOREFRACTION                    = 0.95;      # taxa/genome occupancy for soft-core genes/clusters http://www.biomedcentral.com/1471-2164/13/577/abstract

# paths
our $WORKING_DIR   = '';
our $TMP_DIR       = '';

# files
our $lockfile                            = ".lock";
our $all_fa_file                         = "all.fa";
our $selected_genomes_file               = "selected.genomes"; # records which genomes were previously used to create $bpo_file
our $blast_file                          = "all.blast";
our $bpo_file                            = "all.bpo";
our $blast_file_nr                       = "all.blast.nr";
our $bpo_file_nr                         = "all.bpo.nr";
our $pfam_file                           = "all.pfam";
our $redundant_file                      = "all.redundant";
our $taxa_index_file                     = "taxa_bpo.idx";
our $matrix_file                         = "all_ortho.mtx";
our $mcl_file                            = "all_ortho.mcl";
our $paranoid_file                       = "all_paranoid.out";
our $parameter_BDBH_log                  = "parameter.BDBH.log";
our $parameter_INP_log                   = "parameter.INP.log";
our $parameter_HOM_log                   = "parameter.HOM.log";
our $parameter_OMCL_log                  = "parameter.OMCL.log";
our $parameter_PRND_log                  = "parameter.PRND.log";
our $parameter_COGS_log                  = "parameter.COGS.log";
our $parameter_ISOS_log                  = "parameter.ISOS.log";
our $p2ofilename                         = 'all.p2o.csv'; # files for COG algorithm
our $lsefilename                         = 'all.lsejob.csv';
our $lseoutfilename                      = 'all.lse.csv';
our $coglogfilename                      = 'all.cog.clusters.log';
our $cogblasthits                        = 'hits.csv';

# global variables storing data
our @taxa=();
our %gindex=();     	   # taxon -> [first line, last line, total lines]
our @gindex2=();        # gene_id -> taxon_id
our %blastquery=();
our %pfam_hash;
our %redundant;

our @graph=();
our $last_graph_item = 0;
our %weight=();
our %taxa_bpo_index=(); # hash of arrays[2] that contains first and last file addresses for taxa in $bpo_file

our $blastquery_ref='';

# This subroutine is used to construct the directories to store
# the intermediate files and the final files.
# Arguments: 1 (string) name of desired directory
# Returns:  boolean, 1 if successful, else 0
sub constructDirectory
{
  my ($dirname) = @_;

  $WORKING_DIR = $dirname . '/';
  if(!-e $WORKING_DIR)
  { 
    mkdir($WORKING_DIR) || return 0
  }

  $TMP_DIR = $WORKING_DIR."tmp/";
  if(!-e $TMP_DIR){ mkdir($TMP_DIR); }

  $lockfile                = $TMP_DIR.$lockfile;
  $all_fa_file             = $TMP_DIR.$all_fa_file;
  $selected_genomes_file   = $TMP_DIR.$selected_genomes_file;
  $blast_file              = $TMP_DIR.$blast_file;
  $bpo_file                = $TMP_DIR.$bpo_file;
  $blast_file_nr           = $TMP_DIR.$blast_file_nr;
  $bpo_file_nr             = $TMP_DIR.$bpo_file_nr;
  $pfam_file               = $TMP_DIR.$pfam_file;
  $redundant_file          = $TMP_DIR.$redundant_file;
  $taxa_index_file         = $TMP_DIR.$taxa_index_file;
  $parameter_BDBH_log      = $TMP_DIR.$parameter_BDBH_log;
  $parameter_INP_log       = $TMP_DIR.$parameter_INP_log;
  $parameter_HOM_log       = $TMP_DIR.$parameter_HOM_log;
  $parameter_OMCL_log      = $TMP_DIR.$parameter_OMCL_log;
  $parameter_PRND_log      = $TMP_DIR.$parameter_PRND_log;
  $parameter_COGS_log      = $TMP_DIR.$parameter_COGS_log;
  $parameter_ISOS_log      = $TMP_DIR.$parameter_ISOS_log;
  $matrix_file             = $TMP_DIR.$matrix_file;
  $mcl_file                = $TMP_DIR.$mcl_file;
  $paranoid_file           = $TMP_DIR.$paranoid_file;
  $p2ofilename             = $TMP_DIR.$p2ofilename;
  $lsefilename             = $TMP_DIR.$lsefilename;
  $lseoutfilename          = $TMP_DIR.$lseoutfilename;
  $coglogfilename          = $TMP_DIR.$coglogfilename;
  $cogblasthits            = $TMP_DIR.$cogblasthits;

  return 1;
}

# change default values of $blast_file and $bpo_file as required
sub rename_blast_outfiles
{
  my ($type) = @_;

  if(!$type){ die "# rename_blast_outfiles: need a parameter such as 'features', exit\n"; }
  $blast_file =~ s/all\.blast/all$type\.blast/;
  $bpo_file =~ s/all\.bpo/all$type\.bpo/;
}

# used to put in $parameter_log_file the set of parameters used in a given process, so that they can be checked later
# Two arguments:
# 1. String Variable: type of algorithm for this run
# 2. hash variable: is filled with parameters
# Bruno, sept2011
sub save_params
{
  my($type,%params) = @_;

  my $param_file;

  if($type eq 'BDBH'){ $param_file = $parameter_BDBH_log }
  elsif($type eq 'INP'){ $param_file = $parameter_INP_log }
  elsif($type eq 'HOM'){ $param_file = $parameter_HOM_log }
  elsif($type eq 'OMCL'){ $param_file = $parameter_OMCL_log }
  elsif($type eq 'COGS'){ $param_file = $parameter_COGS_log }
  elsif($type eq 'ISO'){ $param_file = $parameter_ISOS_log }
  else{ $param_file = $parameter_PRND_log }

  open(PARAMS,">$param_file") || die "# save_params : cannot create $param_file\n";

  foreach my $par (sort(keys(%params))){ print PARAMS "$par\t$params{$par}\n"; }

  close(PARAMS);
}

# reads global $parameter_log_file and compares the parameters stored there with those passed in a hash
# returns 1 if params are different
# Two arguments:
# 1. String Variable: type of algorithm for this OMCL|PRND|COGS run
# 2. hash variable: filled with parameters in previous subroutine
# Bruno, sept2010-4
sub check_different_params
{
  my($type,%params) = @_;

  my ($param_file,%old_params);

  if($type eq 'BDBH'){ $param_file = $parameter_BDBH_log }
  elsif($type eq 'INP'){ $param_file = $parameter_INP_log }
  elsif($type eq 'HOM'){ $param_file = $parameter_HOM_log }
  elsif($type eq 'OMCL'){ $param_file = $parameter_OMCL_log }
  elsif($type eq 'COGS'){ $param_file = $parameter_COGS_log }
  elsif($type eq 'ISO'){ $param_file = $parameter_ISOS_log }
  else{ $param_file = $parameter_PRND_log }

  open(OLDPARAMS,$param_file) || die "# check_different_params : cannot read $param_file\n";
  while(<OLDPARAMS>){ if(/^(\S+)\t(\S+)/){ $old_params{$1} = $2 } }
  close(OLDPARAMS);

#if(scalar(keys(%params)) != scalar(keys(%old_params))){ return 1 } # commented oct2011

  foreach my $par (keys(%params)){ if($params{$par} ne $old_params{$par}){ return 1 } }

  return 0;
}

# Three arguments:
# 1. String Variable: Fasta file names, each representing one species, separated by comma, e.g. "Bme.fa,Ccr.fa,Eco.fa"
# 2. String Variable: Name of All_Fasta file
# 3. Boolean: when set if will save sequence length in a returned hash ref
# Nov2011, fills global @taxa, %gindex, @gindex2
# Oct2015, only saves seq lengths if requested; %gindex is reduced to a hash of triplets [1st id, last id, total]
# Nov2015: pre-allocate $n_of_seqs
sub constructAllFasta
{
  my ($fa_files,$n_of_seqs,$all_fa_file,$save_seq_length) = @_;
  my (%seq_len,$id,$len);

  my @fafiles=split (",",$fa_files);
  my $totalgeneno=0;

  $#gindex2 = $n_of_seqs;

  if($all_fa_file){ open (ALLFA,">$all_fa_file") or die "# constructAllFasta: cannot write to $all_fa_file"; }

  foreach my $fafile (sort @fafiles)
  {
    open (FA,$WORKING_DIR.$fafile) or die "# constructAllFasta: cannot open $WORKING_DIR.$fafile";
    $/='>';

    my $taxon=(split (/\.$FASTAEXTENSION$/,$fafile))[0]; # in case original file was .fasta

    push (@taxa,$taxon);
    while (<FA>)
    {
      next unless ($_ ne '>');
      $_=~s/\>$//;
      $_=~s/\r|\n$//;
      if($all_fa_file){ print ALLFA ">$_\n"; } 
      my @lines=split(/\r|\n/,$_);
      $id=shift(@lines);

      if(!defined($gindex{$taxon})){ $gindex{$taxon}[0] = $id }
		  else{ $gindex{$taxon}[1] = $id }
		  $gindex{$taxon}[2]++;	

      next if(!$save_seq_length);
      $len=length(join('',@lines));
      $seq_len{$id}=$len;
    }
    close (FA);
    $/="\n";
    foreach $id ($gindex{$taxon}[0] .. $gindex{$taxon}[1]){ $gindex2[$id]=$taxon }
    $totalgeneno += $gindex{$taxon}[2];
  }

  if($all_fa_file){ close (ALLFA); }

  if($all_fa_file){ print("\n# created file $all_fa_file (".scalar(@fafiles)." genomes, $totalgeneno sequences)\n\n"); }
  else{ print ("\n# ".scalar(@fafiles)." genomes, $totalgeneno sequences\n\n"); }

  return \%seq_len;
} ## constructAllFasta

# One Argument:
# 1. String Variable: fasta file name
# Bruno, Dec2016 (copied & changed named from transcripts.pm)
sub executeDIAMONDMAKEDB 
{
  my ($in) = @_;

  my ($command);

  print("\n# running diamond makedb with $in\n");
  if(-s $in)
  {
    $command = "$DIAMONDMAKEDB --in $in --db $in"; 

    open(EXE,"$command |") || die "# ERROR (executeDIAMONDMAKEDB) : cannot run $command : $!\n";
    while(<EXE>){}
    close(EXE);
  }
  else{ die "# executeDIAMONDMAKEDB : cannot find input FASTA file $in\n"; }
}

# Prepares a string with the right syntax to call DIAMOND blastp  with 1 thread
# Up to 5 arguments:
# 1. String Variable: fasta file name
# 2. String Variable: blast out file name
# 3. String Variable: database name
# 4. String Variable: E-value cutoff
# 5. (Optional) String Variable: hits/query to show in BLAST output
# NOTE: segments masked with --seg yes affect returned %identity
# Bruno, Dec2016
sub format_DIAMONDblastp_command 
{
  my ($infile,$outfile,$db,$Evalue,$hits_to_show) = @_;
  
  my $command = "$DIAMONDPEXE -q $infile --evalue $Evalue -d $db -o $outfile " .
    "--quiet --more-sensitive --seg yes --dbsize $BLAST_DB_SIZE --threads 1 --tmpdir $TMP_DIR "; 

  if($hits_to_show){ $command .= "--max-target-seqs $hits_to_show "; }
  
  $command .= '--outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore ';

  return $command;
}


# One Argument:
# 1. String Variable: fasta file name
# -o t option is removed to solve gene id discrepancy
# Bruno, Jan2011: works with blast and blast+ binaries
sub executeFORMATDB
{
  my ($in,$is_dna,$silent) = @_;

  my ($command);

  print "\n# running makeblastdb with $in\n" if(!$silent);
  if(-s $in)
  {
    if($FORMATDB =~ /makeblastdb/) # blast+
    {
      $command = "$FORMATDB -in $in";
      if($is_dna){ $command .= ' -dbtype nucl ' }
      else{ $command .= ' -dbtype prot ' }
    }
    else
    {
      $command = "$FORMATDB -i $in";
      if($is_dna){ $command .= ' -p F ' }
    }

    open(EXE,"$command |") || die "# ERROR (executeFORMATDB) : cannot run $command : $!\n";
    while(<EXE>){}
    close(EXE);
  }
  else{ die "# executeFORMATDB : cannot find input FASTA file $in\n"; }
} 

# prepares a string with the right syntax to call BLASTP
# Up to 5 params:
# 1. Scalar: fasta file name
# 2. Scalar: blast out file name
# 3. Scalar: database name
# 4. Scalar: pvalue cutoff
# 5. Scalar: max number of hits to report
# Bruno, Feb2016
# Do not mask sequences as we want the alignment
sub format_BLASTP_command_aligns
{
  my ($infile,$outfile,$db,$Evalue,$maxhits) = @_;

  my $command = "$BLASTP -dbsize $BLAST_DB_SIZE " .
      #"-seg yes -soft_masking true " . 
      "-query $infile -evalue $Evalue -db $db -max_target_seqs $maxhits -out $outfile ";

  return $command;
}

# Same params as previous sub, plus 
# 6. Scalar: BLASTN algorithm (optional)
# Do not mask sequences as we want the alignment
sub format_BLASTN_command_aligns
{
  my ($infile,$outfile,$db,$Evalue,$maxhits,$task) = @_;

  my $command = "$BLASTN -dbsize $BLAST_DB_SIZE " . 
      #"-soft_masking true " . #"-outfmt 4 " .
      "-query $infile -evalue $Evalue -db $db -max_target_seqs $maxhits -max_hsps 1 -out $outfile ";

  if($task){ $command .= "-task $task "; }
  else{ $command .= "-task megablast "; }

  return $command;
}

# prepares a string with the right syntax to call BLASTP
# Up to 5 arguments:
# 1. String Variable: fasta file name
# 2. String Variable: blast out file name
# 3. String Variable: database name
# 4. String Variable: pvalue cutoff
# 5. (Optional) String Variable: hits/query to show in BLAST output
# 6. (Optional) String: run mode
# Bruno, Oct2015
sub format_BLASTP_command
{
  my ($infile,$outfile,$db,$Evalue,$hits_to_show,$scapequotes) = @_;

  #if($BLASTP =~ /blastall/){
  #  #$command = "$BLASTP -p blastp -a $BLAST_NOCPU -z $BLAST_DB_SIZE -m 8 ".
  #  $command = "$BLASTP -p blastp -z $BLAST_DB_SIZE -m 8 ". # standard settings
  #    "-i $infile -e $Evalue -d $db -o $outfile ";
  #  if($scapequotes){ $command .= "-F \\'m\\ S\\' " } # filter as in PMID:18042555
  #  else{ $command .= "-F 'm S' " }
  #  if($hits_to_show){ $command .= "-v $hits_to_show -b $hits_to_show"; } }

  my $command = "$BLASTP -dbsize $BLAST_DB_SIZE " .
      "-seg yes -soft_masking true " . # filter as in PMID:18042555
      "-query $infile -evalue $Evalue -db $db -out $outfile ";

  # COGs expects 12-column output, but actually only parses 0,1,6-11, so we added qlen and slen as 4 & 5
  if($scapequotes)
  {
    $command .= " -outfmt \\'6\\ qseqid\\ sseqid\\ pident\\ length\\ qlen\\ slen\\ qstart\\ qend\\ sstart\\ send\\ evalue\\ bitscore\\' "
  }
  else
  {
    $command .= " -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' "
  }

  if($hits_to_show){ $command .= "-max_target_seqs $hits_to_show"; }

  return $command;
}

# Bruno, Oct2015
sub format_BLASTN_command
{
  my ($infile,$outfile,$db,$Evalue,$hits_to_show,$scapequotes,$task) = @_;

  #if($BLASTN =~ /blastall/){
  #  #$command = "$BLASTN -p blastn -a $BLAST_NOCPU -z $BLAST_DB_SIZE -m 8 " .
  #  $command = "$BLASTN -p blastn -z $BLAST_DB_SIZE -m 8 " . # standard settings
  #    "-i $infile -e $Evalue -d $db -o $outfile ";
  #  if($hits_to_show){ $command .= "-v $hits_to_show -b $hits_to_show"; }}

  my $command = "$BLASTN -dbsize $BLAST_DB_SIZE " . 
      "-soft_masking true " .
      "-query $infile -evalue $Evalue -db $db -out $outfile ";

  if($scapequotes)
  { 
    $command .= " -outfmt \\'6\\ qseqid\\ sseqid\\ pident\\ length\\ qlen\\ slen\\ qstart\\ qend\\ sstart\\ send\\ evalue\\ bitscore\\' " 
  } 
  else
  { 
    $command .= " -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' " 
  }

  if($task){ $command .= "-task $task "; }
  else{ $command .= "-task megablast "; }
  if($hits_to_show){ $command .= "-max_target_seqs $hits_to_show "; }

  return $command;
}

# Bruno Mar2013
sub format_SPLITBLAST_command
{
  return $ENV{"EXE_SPLITBLAST"}." $BLAST_NOCPU ";
}

# Bruno Mar2013
sub format_HMMPFAM_command
{

  #return $ENV{"EXE_HMMPFAM"}." --cpu $BLAST_NOCPU ".$ENV{"PFAMDB"}.' ';
  return $ENV{"EXE_HMMPFAM"}.' --cpu 1 '.$ENV{"PFAMDB"}.' ';
}

sub format_SPLITHMMPFAM_command
{
  return $ENV{"EXE_SPLITHMMPFAM"}." $BLAST_NOCPU ";
}

# Up to four Arguments:
# 1. String Variable: matrix file name
# 2. String Variable: mcl file name
# 3. Number Variable: inflation parameter for MCL algorithm
# 4. (Optional) String Variable: verbose flag
# Modified by Bruno Dec2011
sub executeMCL
{
  my ($matrix_file, $mcl_file, $inflation, $verbose)  = @_;

  if($verbose){ $verbose = '' }
  else{ $verbose = " > /dev/null" }

  print "\n# running MCL (inflation=$inflation) ...\n";

  #die "$MCL $matrix_file -I $inflation -o $mcl_file --silent $verbose 2>&1 "; # debug
  #system("$MCL $matrix_file -I $inflation -o $mcl_file --silent $verbose 2>&1 ");
  system("$MCL $matrix_file -I $inflation -o $mcl_file -V all -te $BLAST_NOCPU $verbose 2>&1 ");
  if($? != 0)
  {
    die "# ERROR: failed running MCL (probably ran out of memory, option -s recommended)\n";
  }
  elsif(!-s $mcl_file)
  {
    die "# ERROR: failed generating $mcl_file file (OMCL)\n";
  }
  else
  {
    print("# running MCL finished\n");
  }

}

# save to file hash of redundant isoforms (global %redundant) for later queries from check_BDBHs.pl
# Bruno Jan2015
sub write_redundant_hash
{
  my ($outfile) = @_;

  open(REDFILE,">$outfile") || die "# write_redundant_hash : cannot create $outfile $!\n";

  foreach my $isof (keys(%redundant))
  {
    print REDFILE "$isof\t$redundant{$isof}\n";
  }

  close(REDFILE);
}

#Query sequence: 11
#[..]
#Parsed for domains:
#Model      Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
#--------   ------- ----- -----    ----- -----      -----  -------
#PF08534.2    1/1       3   161 ..     1   179 []    29.5  1.4e-05
#PF00578.13   1/1       4   140 ..     1   132 []   163.5  6.2e-46
#PF10417.1    1/1     150   193 ..     1    63 []    78.2  3.1e-20
#//
# Added by Bruno March2010, works with HMMER2.3&3
# Updated Mar2016
sub pfam_parse
{
  my ( $parseoutfile, @pfam_files ) = @_;

  my ($n_of_lines,$id,$family,$from,$to,$Evalue) = (0);

  open(PARSEDPFAM,">$parseoutfile") || die "# pfam_parse : cannot create $parseoutfile $!\n";

  foreach my $pfamfile (@pfam_files)
  {
    my (%data,%descr);
    open(PFAM,$pfamfile) || die "# pfam_parse : cannot read $pfamfile\n"; #print "# $pfamfile\n";
    while(<PFAM>)
    {
      if(/^Query sequence:\s+(\d+)/ || /^Query:\s+(\d+)/){ $id = $1; } # works with HMMER2.3&3
      elsif(/^(PF\S+)\.\d+\s+\S+\s+(\d+)\s+(\d+)\s+\S+\s+\d+\s+\d+\s+\S+\s+\S+\s+(\S+)/) # works with HMMER2.3
      {
        # read all PFAM domains reported and capture Evalue and coordinates
        ($family,$from,$to,$Evalue) = ($1,$2,$3,$4);
        $data{$id}{$Evalue} = "$family,$from,$to"; #print "$id $family $from $to $Evalue\n";
      }
      elsif(/^>> (PF\S+)\.\d+\s+(.*?)\n/) # works with HMMER3
      { 
        #>> PF10602.5  26S proteasome subunit RPN7
        $family = $1; 
        $descr{$family} = $2; 
      }  
      elsif(/^\s+\d+\s+[!|\?]\s+/)# works with HMMER3
      {
        #>> PF00132.17  Bacterial transferase hexapeptide (three repeats) <== OJO, FUSIONA REPEATS
        #   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
        #   1 !   15.7   0.4   2.6e-06    0.0078       1      16 [.     283     298 ..     283     300 .. 0.93
        #   1 !   98.5   0.1     4e-32   2.4e-28      # hmmscan :: search sequence(s) against a profile database
        my @data = split;
        ($Evalue,$from,$to) = ($data[5],$data[9],$data[10]);
        next if($from eq 'search');
        $data{$id}{$Evalue} = "$family,$from,$to"; #print "$id $family $from $to $Evalue\n";
      }
      elsif(/no hits above thresholds/ || /No hits detected that satisfy reporting thresholds/)# works with HMMER2.3&3
      {
        $data{$id} = 'nohits';
      }
    }
    close(PFAM);

    # for each annotated protein sort PFAM domains in terms of Evalue and take best alignment
    # (in order to cope with cases where two PFAM domains cover roughly the same sequence)
    foreach $id (sort {$a<=>$b} keys(%data))
    {
      my ($dom_string,$descr_string,%PFAMcover,$midpoint,$domain,$overlapOK) = ('');

      if($data{$id} eq 'nohits')
      {
        print PARSEDPFAM "$id\t\t\n";
        next;
      }

      foreach $Evalue (sort {$a<=>$b} keys(%{$data{$id}}))
      {
        $overlapOK = 0;
        ($family,$from,$to) = split(/,/,$data{$id}{$Evalue});
        $midpoint = $from+(($to-$from)/2); #print "## $id $Evalue $from $to $midpoint $family\n";
        foreach $domain (keys(%PFAMcover))
        {
          if($midpoint > $PFAMcover{$domain}{'from'} && $midpoint < $PFAMcover{$domain}{'to'})
          {
            $overlapOK = 1;
          }
        }
        next if($overlapOK);
        $PFAMcover{$family.'_'.$midpoint}{'from'} = $from;
        $PFAMcover{$family.'_'.$midpoint}{'to'} = $to; # Evalue is not need anymore
      }

      foreach $domain (sort {$PFAMcover{$a}{'from'}<=>$PFAMcover{$b}{'from'}} keys(%PFAMcover))
      {
        $family = (split(/_/,$domain))[0];

        # print ">$id $family $PFAMcover{$domain}{'from'} $PFAMcover{$domain}{'to'}\n";
        $dom_string   .= "$family,";
        $descr_string .= "$descr{$family};" || 'NA;';
      }

      print PARSEDPFAM "$id\t$dom_string\t$descr_string\n";
    }
  }

  close(PARSEDPFAM);

  return $n_of_lines;
}

# sorts and merges partial BLAST/DIAMOND results using GNU sort
# optionally also compresses results as a side-effect
# arguments:
# 1. String Variable: blast out file
# 2. boolean flag stating whether secondary hsps are to be conserved 
# 3. boolean flag stating whether BLAST outfiles are to be compressed
# 4. array of sorted blast output filenames
# uses globals: $TMP_DIR , $SORTLIMITRAM, $SORTBIN, $GZIPBIN
# Updated Mar2019
sub sort_blast_results
{
  my ($sorted_outfile,$keep_secondary_hsps,$compress_blast,@blast_outfiles) = @_;

  my ($file,$size,$files,$n_of_files,$root,$cleantmpfile,$sortedtmpfile);
  my (@tmpfiles,@roots,%cluster);

  # check max file descriptors in OS
  my $max_os_files = qx{echo `ulimit -n`};

  if(-s $sorted_outfile){ unlink($sorted_outfile) }

  foreach $file (@blast_outfiles)
  {
    # _homologues/_Buch_aph_APS.faa.fasta_Buch_aph_APS.faa.fasta.blast
    $root = (split(/\//,$file))[-1];
    $cleantmpfile = $TMP_DIR . "/_tmp_$root";
    $root = (split(/\.fasta_/,$root))[0];
    if(!grep(/$root/,@roots)){ push(@roots,$root); } # remember roots in proper order
    if($file =~ /$root\.fasta$root/)
    {
      # make sure inparalogs are up in the list
      unshift(@{$cluster{$root}},$cleantmpfile);
    }
    else{ push(@{$cluster{$root}},$cleantmpfile); }
    
    if($compress_blast)
    {
      if(!-s $file.'.gz')
      {
        system("$GZIPBIN $file");
        if($? != 0)
        {
          die "# sort_blast_results : cannot compress $file\n";
        }
      }
      
      open(ORIG,"$GZIPBIN -dc $file |") || die "# sort_blast_results : cannot read $file.gz\n";
    }
    else
    {
      open(ORIG,$file) || die "# sort_blast_results : cannot read $file\n";
    }

    # conserve secondary hsps attached to first one to avoid breaking ordered BLAST output
    # 90695	5112	70.57	1043	176	17	1	940	248	1262	0.0	1358
    # 90695	5112	24.88	422	223	16	183	544	229	616	2e-09	65.5
    # 90695	5279	60.49	1030	286	16	1	939	248	1247	0.0	1136
    my ($pqry,$psbj,$qry,$sbj) = ('-1','-1');
    open(CLEANBLAST,">$cleantmpfile") || die "# sort_blast_results : cannot create $cleantmpfile\n";
    
    while(<ORIG>)
    {
      chomp;
      if(/^(\d+)\t(\d+)/){ ($qry,$sbj) = ($1,$2) }
      
      if($keep_secondary_hsps) # Mar2015 conserves secondary hsps 
      {
        if($qry ne $pqry || $sbj ne $psbj){ print CLEANBLAST "\n$_"; } # adds empty line on top of file
        else{ print CLEANBLAST "\thsp\t$_"; } # concat secondary hsps on same line 
      }
      else # original version, only keeps primary hsp
      {
        if($qry ne $pqry || $sbj ne $psbj){ print CLEANBLAST "$_\n"; } 
      }  
      
      ($pqry,$psbj) = ($qry,$sbj);
    }
    close(ORIG);
    
    if($keep_secondary_hsps){ print CLEANBLAST "\n"; } # add last newline if required
    
    close(CLEANBLAST);
   
    push(@tmpfiles,$cleantmpfile);      
  }

  foreach $root (@roots)
  {
    ($files,$size,$n_of_files) = ('',0,0);
    foreach $file (@{$cluster{$root}})
    {
      $size += (-s $file) / (1024*1024);
      $files .= " $file";
		$n_of_files++;
    }
    printf("# sorting %s results (%1.2gMB)\n",$root,$size);

    # deprecated as it re-orders what blast had already ordered and requires more temporary files July2010
    #$tmpfile = $TMP_DIR . "/_tmp_$root.blast";
    #foreach $file (@{$cluster{$root}}){ system("cat $file >> $tmpfile"); }
    #system("$SORTBIN -T=$TMP_DIR -s -g -k1 -k11 $tmpfile > $sortedtmpfile"); # order by queryID > E-value
    $sortedtmpfile = $TMP_DIR . "/_tmp_$root.blast.sorted";

    # http://www.gnu.org/software/coreutils/faq/coreutils-faq.html #Sort-does-not-sort-in-normal-order_0021
    $ENV{'LC_ALL'} = 'POSIX';
   
    #hits are sorted by query ID and E-value
    #hits with same E-value and diff bitscore will be taken care later 
    #asumes BLAST results are ordered
    #tries to read & merge all individual files at once to avoid temp files
    #https://stackoverflow.com/questions/6598573/how-to-merge-sorted-files-without-using-a-temporary-file
    if($n_of_files > $max_os_files){ $n_of_files = $max_os_files }
    system("$SORTBIN --temporary-directory=$TMP_DIR --batch-size=$n_of_files -s -k1g -k11g -m $files > $sortedtmpfile");
    if(!-s $sortedtmpfile) #shell sort failed
    {
      merge_BLAST_files($sortedtmpfile,@{$cluster{$root}});
      if(!-s $sortedtmpfile)
      {
        die "# sort_blast_results : failed while generating partial $sortedtmpfile\n";
      }
    }
    
    if($keep_secondary_hsps)
    {
      open(SORTED,">>$sorted_outfile") || die "# sort_blast_results : cannot open $sorted_outfile in append mode\n"; 
      
      open(TMPSORTED,$sortedtmpfile);
      while(<TMPSORTED>)
      {
        next if(/^$/);
        $_ =~ s/\thsp\t/\n/g;
        print SORTED $_;
      }
      close(TMPSORTED);
      
      close(SORTED);
    }
    else{ system("cat $sortedtmpfile >> $sorted_outfile") }

    unlink($sortedtmpfile);
  }
  if(!-s $sorted_outfile){ die "# sort_blast_results : failed while generating $sorted_outfile\n"; }
  unlink(@tmpfiles);
}

sub merge_BLAST_files
{
  # Adapted from File-Sort-1.01(http://search.cpan.org/perldoc?File::Sort)
  # Assumes infiles are BLAST output files in tabular format with sequences
  # identified by natural numbers, such as 11,12,1439,1440 in the sample:
  #11 1439    78.24   625 136 0   4   628 6   630 0.0  993
  #12 1440    80.88   272 52  0   1   272 1   272 7e-125   446
  # Order of infiles is RELEVANT as merging is stable, so sequences from
  # the first files will be given sorting priority

  my ($outfile,@infiles) = @_;

  my (%order,@fhorder,$filein,$id,$first,$n_of_fhs,$curr,$line);
  our %fh;

  # for Schwartzian transform (ST) see
  # http://www.hidemail.de/blog/perl_tutor.shtml#sort_orcish_schwartzian
  sub blastsort { $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] }
  sub blastmap
  {
    my @tmp = split(/\s+/,$fh{$_});
    [$_,$fh{$_},$tmp[0],$tmp[10]]

    # returns anonymous array[4]: filehandle, line,query_id,E-value
  }

  ## open all input BLAST files and keep filehandles, ordered 0..N
  $n_of_fhs = 0;
  foreach $filein (@infiles)
  {
    $id = gensym(); # get valid id for filehandle

    open($id,$filein) ||
      die "# merge_BLAST_files : cannot read $filein: $!";

    $order{$n_of_fhs} = $id;
    $n_of_fhs++;
  }
  @fhorder = (0 .. $n_of_fhs-1);

  ## open outfile and ensure IO buffer is used
  $| = 0;
  unlink($outfile) if(-s $outfile);
  open(OUT,">$outfile") ||
    die "# merge_BLAST_files : cannot create $outfile: $!";

  ## get first BLAST line from all filehandles
  %fh = map {
    my $fh = $order{$_};
    ($_ => scalar <$fh>);
    } @fhorder;

  ## start merging BLAST lines
  while(scalar(@fhorder)>1)
  {
    ($first) = (map {$_->[0]} sort blastsort map &blastmap, @fhorder); #ST

    print OUT $fh{$first};

    $curr = $order{$first};
    $line = scalar <$curr>;
    if(defined($line)) # update current filehandle
    {
      $fh{$first} = $line;
    }
    else # exhausted filehandle
    {
      @fhorder = grep { $_ ne $first } @fhorder;
    }
  }

  ## take care of last filehandle left and close file
  print OUT $fh{$fhorder[0]};
  $curr = $order{$fhorder[0]};
  while(<$curr>){ print OUT }
  close(OUT);

  ## close all input files
  foreach $id (0 .. $n_of_fhs-1){ close($order{$id}); }

  if(!-s $outfile){ return 0 }
  else{ return 1 }
}

# One argument:
# 1. String Variable: blast out file
sub blast_parse_COG
{
  my ($blastfile,$skipsort) = @_;

  print("\n# parsing file (COG)\n");
  if(!$skipsort) # if called before sorting blast output, never used in get_homologues.pl, July2010
  {
    print("\n# sorting file (COG)\n");

    # http://www.gnu.org/software/coreutils/faq/coreutils-faq.html
    # #Sort-does-not-sort-in-normal-order_0021
    $ENV{'LC_ALL'} = 'POSIX';
    my $sortedblastfile  = $blastfile . ".sorted";
    system("$SORTBIN -T=$TMP_DIR -k1g -k11g $blastfile > $sortedblastfile"); # order by queryID > E-value
    system("mv -f $sortedblastfile $blastfile"); # save disk space
    print "# done sorting file (COG)\n";
  }

  # clean previous links (crean problemas alternando con -a, Jul2012)
  if(-l $TMP_DIR.'all.blast.tab')
  {
    unlink($TMP_DIR.'all.blast.tab') ||
      die "# ERROR: failed removing symbolic link $TMP_DIR\.all.blast.tab (COG)\n";
  }
  if(-l $TMP_DIR.'allfeatures.blast.tab')
  {
    unlink($TMP_DIR.'allfeatures.blast.tab') ||
      die "# ERROR: failed removing symbolic link $TMP_DIR\.allfeatures.blast.tab (COG)\n";
  }

  # produces hash.csv
  system("$MAKEHASHEXE -i=$p2ofilename -o=$TMP_DIR");
  if($? != 0)
  {
    die "# ERROR: failed producing blast hash (COG)\n"
  }
  elsif(!-e $blastfile.'.tab')
  {

    #$READBLASTEXE requires a .tab extension
    if(!symlink( $blastfile ,  $blastfile.'.tab' ))
    {
      die "# ERROR: failed making symbolic link (COG)\n"
    }
  }

  # produces hits.csv, self.csv, COGreadblast.log, query2subject.csv
  system("$READBLASTEXE -d=$TMP_DIR -u=$TMP_DIR -f=$TMP_DIR -s=$TMP_DIR -q=2 -t=2");
  if($? != 0)
  {
    die "# ERROR: failed parsing file (COG): $?\n"
  }
  elsif(!-e $cogblasthits)
  {
    die "# ERROR: failed producing hits file (COG)\n"
  }
  else
  {
    print("# parsing file (COG) finished\n");
  }
}

# Three Arguments:
# 1. String Variable: blast out file
# 2. String Variable: blast parse out file
# 3. (optional) Flag: skip sort step if requested
# 4. (optional) Flag: trim BLAST multi-hsp overlaps, else over-estimate coverage occasionally
# Q=query,S=subject, assumes Qids are integers
# Updated May2016
sub blast_parse
{
  my ($blastfile,$parseoutfile,$skipsort,$trim) = @_;

  my $n_of_bpo_lines = 0;

  #1	1	100.00	521	521	521	1	521	1	521	0.0	1025
  #1	1058	60.24	503	198	521	515	515	14	515	2e-173	 610
  if(!$skipsort)
  {
    print("\n# sorting file\n");

    # http://www.gnu.org/software/coreutils/faq/coreutils-faq.html
    # #Sort-does-not-sort-in-normal-order_0021
    $ENV{'LC_ALL'} = 'POSIX';
    my $sortedblastfile  = $blastfile . ".sorted";
    system("$SORTBIN -T=$TMP_DIR -g -k1 -k11 $blastfile > $sortedblastfile"); # order by queryID > E-value
    system("mv -f $sortedblastfile $blastfile"); # save disk space
    print("# done sorting file\n");
  }

  my $blastfile_size = sprintf("%1.2gMB",(-s $blastfile) / (1024*1024));

  open (PARSEOUT,">$parseoutfile") || die "# blast_parse : cannot create $parseoutfile\n";
  print("\n# parsing blast result! ($blastfile , $blastfile_size)\n");

  open(BLASTOUT,$blastfile) || die "# blast_parse : cannot read $blastfile\n";

  my ($Qid,$Sid,$percID,$Qaln_len,$Saln_len,$Qlength,$Slength,$Qstart,$Qend,$Sstart,$Send,$Evalue,$bits);
  my ($pQid,$pSid,$ppercID,$pQaln_len,$pSaln_len,$pQlength,$pSlength,$pQstart,$pQend,$pSstart,$pSend,$pEvalue,$pbits); # p=previous
  my ($Qcov,$Scov,$n_of_hsp,$simspan,$aln_length,$midpointQ,$midpointS,$h,$redhsp,@hsps);

  $pQid=$pSid=-1;
  ($n_of_hsp,$simspan,$redhsp) = (1,'',0);

  while(<BLASTOUT>)
  {
    chomp;
    ($Qid,$Sid,$percID,$aln_length,$Qlength,$Slength,$Qstart,$Qend,$Sstart,$Send,$Evalue,$bits) = split(/\t/);
    if($Send<$Sstart){ ($Sstart,$Send) = ($Send,$Sstart) } #print "$Qstart $Qend\n";

    #$Evalue = numeric_pvalue($Evalue); 
    #next if($pv_cutoff && $Evalue > $pv_cutoff); # deprecated
    #$percID = sprintf("%1.0f",$percID); 

    if($Qid eq $pQid && $Sid eq $pSid) # multiple hsp match
    {
		# skip last parsed hsp if redundant 
		$midpointQ = $Qstart + (($Qend-$Qstart)/2);  
		$midpointS = $Sstart + (($Send-$Sstart)/2);
		for($h=0;$h<$n_of_hsp;$h++)
		{ 
      #make sure redundant hsps are ignored and otherwise trim ovelapping ends 
      #to avoid over-estimating cover (bitscore and $ID might still be over-estimated)
		  if( ($midpointQ >= $hsps[$h][0] && $midpointQ <= $hsps[$h][1]) 
			 || ($midpointS >= $hsps[$h][2] && $midpointS <= $hsps[$h][3])
		  ){ $redhsp = 1; last }

		  if($trim)
      {
          #trim overlapping Q coords
          if($Qstart < $hsps[$h][1] && $Qend > $hsps[$h][1]){ $Qstart = $hsps[$h][1]+1 }
          elsif($Qend > $hsps[$h][0] && $Qstart < $hsps[$h][0]){ $Qend = $hsps[$h][0]-1 }
          #trim overlapping S coords
          if($Sstart < $hsps[$h][3] && $Send > $hsps[$h][3]){ $Sstart = $hsps[$h][3]+1 }
          elsif($Send > $hsps[$h][2] && $Sstart < $hsps[$h][2]){ $Send = $hsps[$h][2]-1 }
		  }
		}
		if($redhsp || $Qend<=$Qstart || $Send<=$Sstart)
		{
		  $redhsp = 0;
        next;
		}
		else
		{        
        $Qaln_len = ($Qend-$Qstart)+1;
        $Saln_len = ($Send-$Sstart)+1;
		}

      # update scores after adding this valid hsp
      if($Evalue > $pEvalue){ $Evalue = $pEvalue } # toma mejor evalue como representativo 
      $bits += $pbits;
      $percID += $ppercID;
      $Qaln_len += $pQaln_len;
		  $Saln_len += $pSaln_len;
      push(@hsps,[$Qstart,$Qend,$Sstart,$Send]);
      $simspan .= "$n_of_hsp:$pQstart-$pQend:$pSstart-$pSend.";
      $n_of_hsp++;
    }
    elsif($pQid != -1) # normal,single hsp match, handle previous line
    {
      if($n_of_hsp>1){ $ppercID = sprintf("%1.2f",$ppercID/$n_of_hsp) }

      $simspan .= "$n_of_hsp:$pQstart-$pQend:$pSstart-$pSend.";
		  $Qcov = sprintf("%1.0f",100 * $pQaln_len / $pQlength);
		  $Scov = sprintf("%1.0f",100 * $pSaln_len / $pSlength);

		  # move length, span, bits to the end as they are rarely needed
      print PARSEOUT "$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits\n";
      $n_of_bpo_lines++;
      
		  # initialize
		  ($n_of_hsp,$simspan,$redhsp) = (1,'',0);
      @hsps = ([$Qstart,$Qend,$Sstart,$Send]);
      $Qaln_len = ($Qend-$Qstart)+1;
      $Saln_len = ($Send-$Sstart)+1;
    }
    else # $pQid = -1, first bpo line ever
    {
      @hsps = ([$Qstart,$Qend,$Sstart,$Send]);    
      $Qaln_len = ($Qend-$Qstart)+1;
      $Saln_len = ($Send-$Sstart)+1;
    }

    # store this BLAST lines to be processed after the next one is parsed
    ($pQid,$pSid,$ppercID,$pQaln_len,$pSaln_len,$pQlength,$pSlength,$pQstart,$pQend,$pSstart,$pSend,$pEvalue,$pbits) =
      ($Qid,$Sid,$percID,$Qaln_len,$Saln_len,$Qlength,$Slength,$Qstart,$Qend,$Sstart,$Send,$Evalue,$bits);
  }

  # handle last output line
  if(!$redhsp)
  {
    if($n_of_hsp>1){ $ppercID = sprintf("%1.2f",$ppercID/$n_of_hsp) }
  
    $simspan .= "$n_of_hsp:$pQstart-$pQend:$pSstart-$pSend.";
    $Qcov = int(100 * $pQaln_len / $pQlength);
    $Scov = int(100 * $pSaln_len / $pSlength);
    print PARSEOUT "$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits\n";
    $n_of_bpo_lines++;
  }

  close(BLASTOUT);

  print("# parsing file finished\n");
  close(PARSEOUT);

  return $n_of_bpo_lines;
} ## blast_parse

# Filters results for redundant sequences among a large blast report.
# Redundant sequences are labelled in $ref_rd_sequences
# Bruno Jan2015
sub nr_blast_report
{
  my ($blastfile,$nrblastfile,$ref_rd_sequences) = @_;

  my ($nr_lines) = (0);

  open (NRBLAST,">$nrblastfile") || die "# nr_blast_parse : cannot create $nrblastfile\n";

  open(BLASTOUT,$blastfile) || die "# nr_blast_parse : cannot read $blastfile\n";

  print("\n# redundancy-filtering blast file\n");

  while(<BLASTOUT>)
  {
    if(/^(\d+)\t(\d+)/)
    {
      next if($ref_rd_sequences->{$1} || $ref_rd_sequences->{$2});
      print NRBLAST;
      $nr_lines++;
    }
  }

  close(BLASTOUT);
  close(NRBLAST);

  print("# created nr blast file\n");

  return $nr_lines;
}

# This module is used to open BPO (Blast Parse Out) file,
# to provide an filehandle for later use in blast query process.
# One Argument:
# 1. String Variable: BPO file name
# Last modified: 07/21/04
sub open_bpofile
{
  my $bpo_file = $_[0];
  open (BPOFILE,$bpo_file) || die("# open_bpofile : can't open $bpo_file file");
} # open_bpofile

# Make pvalue numeric, used by subroutine blast_parse
# One Arguments:
# 1. String Variable: pvalue
# Last modified: 07/19/04
#sub numeric_pvalue {
#  my $p=$_[0];
#  if ($p=~/^e-(\d+)/){ return '1e-'.$1; }
#  else {return $p}
#} # numeric_pvalue

# fills global %redundant with redundancies read from $infile, which usually is $redundant_file
# Updated by Bruno Oct2015
sub construct_redundant_hash
{
  my ($infile,$target,$verbose) = @_;
  open(RED,$infile) || die "# construct_redundant_hash : cannot read $infile : $!\n";
  while(<RED>)
  {
    if(/^(\S+)\t(\S+?)\n/)
    {
      next if(@gindex2 && !defined($gindex2[$1])); # read only sequences that belong to valid taxa, and save memory
      next if($target && $target ne $1);
      $redundant{$1} = $2; #print "$1\t$2\n";
    }
  }
  close(RED);

  print "# construct_redundant_hash: records = ".scalar(keys(%redundant))."\n" if($verbose);
}

# fills global %pfam_hash with Pfam assignments parsed from $infile, which usually is $pfam_file
# Updated Bruno Feb2016
sub construct_Pfam_hash
{
  my ($infile,$verbose) = @_;
  open(PARSEDPFAM,$infile) || die "# construct_Pfam_hash : cannot read $infile : $!\n";
  while(<PARSEDPFAM>)
  {
    if(/^(\S+)\t(\S+?)\t/)
    {
      next if(@gindex2 && !defined($gindex2[$1])); # read only sequences that belong to valid taxa, and save memory
      $pfam_hash{$1} = $2; #print "$1\t$2\n";
    }
  }
  close(PARSEDPFAM);

  print "# construct_Pfam_hash: records = ".scalar(keys(%pfam_hash))."\n" if($verbose);
}

# Reads $infile, which usually is $pfam_file, and counts the occurrence of Pfam domains in sequences
# indicated in arrays $ref_ids1 and $ref_ids2. Usually 1 will be a subset of 2, but it is not required.
# Only the fist occurrence of a domain in a sequence is counted.
# If optional $ref_clusters is passed, domain occurrences are counted once per cluster
# Bruno Mar-Sept2016
sub parse_Pfam_freqs
{
  my ($infile,$ref_ids1,$ref_ids2,$ref_clusters) = @_;
  
  my ($pfam,$domains,$id,$descriptions,$cluster);
  my (%counts1,%counts2,%full_text,%cluster_pfam);
  my ($idx1,$idx2) = (0,0);
  
  #printf("# parse_Pfam_freqs: last1 = %d last2 = %d sequence ids\n\n",
  #  $ref_ids1->[$#{$ref_ids1}],$ref_ids2->[$#{$ref_ids2}]);
  
  open(PARSEDPFAM,$infile) || die "# parse_Pfam_freqs : cannot read $infile : $!\n";
  while(<PARSEDPFAM>)
  {
    #31   
    #32	PF03602,	Conserved hypothetical protein 95;
    #33	PF02881,PF00448,	SRP54-type protein, helical bundle domain;SRP54-type protein, GTPase domain;
    chomp;
    ($id,$domains,$descriptions) = split(/\t/,$_); 
    
    last if($id > $ref_ids1->[$#{$ref_ids1}] && $id > $ref_ids2->[$#{$ref_ids2}]);
    
    # update lagging indexes
    while($ref_ids1->[$idx1] < $id && $idx1 < $#{$ref_ids1}){ $idx1++ }
    while($ref_ids2->[$idx2] < $id && $idx2 < $#{$ref_ids2}){ $idx2++ }
    #print "id:$id $ref_ids1->[$idx1] $ref_ids2->[$idx2] $descriptions\n";

    # skip ids without Pfam annotations     
    next if($domains eq '');
    
    # check cluster to which this sequence belongs
    if($ref_clusters)
    { 
      $cluster = $ref_clusters->{$id}; #print "$id -> $cluster\n";
    }
    
    # split domain and description strings
    my %seen;
    my @doms   = split(/,/,$domains);
    my @descrs = split(/;/,$descriptions);
    
    # make non-redundant hash of Pfam domains and save Pfam text descriptions
    foreach $pfam (0 .. $#doms)
    {
      next if($seen{$doms[$pfam]});
      
      # if required make sure that domains are counted once per cluster
      next if($ref_clusters && defined($cluster_pfam{$ref_clusters->{$id}}{$doms[$pfam]}));
      
      $seen{$doms[$pfam]}++;
      
      if(!$full_text{$doms[$pfam]})
      {
        $full_text{$doms[$pfam]} = $descrs[$pfam];
        #print ">$doms[$pfam]<>$descrs[$pfam]<\n";
      }
      
      if($ref_clusters)
      { 
        # add this Pfam domain to cluster annotation
        $cluster_pfam{$ref_clusters->{$id}}{$doms[$pfam]}=1; 
      }
    }
    
    # add individual domains
    if($id == $ref_ids1->[$idx1])
    {
      foreach $pfam (keys(%seen)){ $counts1{$pfam}++ }
      if($idx1 < $#{$ref_ids1}){ $idx1++ }
    }
    
    # add individual domains
    if($id == $ref_ids2->[$idx2])
    { 
      foreach $pfam (keys(%seen)){ $counts2{$pfam}++ }  
      if($idx2 < $#{$ref_ids2}){ $idx2++ }
    } #printf("# %d %d\n",scalar(keys(%counts1)),scalar(keys(%counts2)));
  }
  close(PARSEDPFAM);

  printf("# parse_Pfam_freqs: set1 = %d Pfams set2 = %d Pfams\n\n",
    scalar(keys(%counts1)),scalar(keys(%counts2)));

  # fill blanks %counts1 to obtain equal tables
  foreach $pfam (keys(%counts2))
  {
    next if(defined($counts1{$pfam}));
    $counts1{$pfam}=0;
  }
  
  # detect errors
  #foreach $pfam (sort {$counts1{$b}<=>$counts1{$a}} keys(%counts1)){ print "$pfam $counts1{$pfam} $counts2{$pfam}\n" }
  
  return (\%counts1,\%counts2,\%full_text);
}




# finds file coordinates for taxa as blast queries and saves them in global %taxa_bpo_index
# Updated Bruno Oct2015
sub construct_taxa_indexes
{
  my ($bpo_file) = @_;

  my ($n_of_taxa,$add,$taxa,%taxa_add)=(0);

  if(!%gindex) # try to re-use previously calculated index (in master node)
  {
    if(!-s $taxa_index_file){ die("# construct_taxa_indexes: cannot find $taxa_index_file\n") }
    my ($ref_data_structure); #print "recovering $taxa_index_file\n";
    eval { $ref_data_structure = retrieve($taxa_index_file) };
    if(!$@) #  compatible previous binary (avoid trans-arch errors)
    {
      if(defined($ref_data_structure))
      {
        %taxa_bpo_index = %$ref_data_structure;
      }
      else{ die("# construct_taxa_indexes: $taxa_index_file is empty\n") }
    }
    else{ die("# construct_taxa_indexes: cannot recover $taxa_index_file : $@\n") }
  }
  else
  {
    open(BPOFILE,$bpo_file)|| die("# construct_taxa_indexes: cannot open $bpo_file file");
    $add = tell(BPOFILE);
    while (<BPOFILE>)
    {
      #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      if(/^(\d+)/)
      {
        $taxa = $gindex2[$1];
        if(!defined($taxa_add{$taxa}[0])){ $taxa_add{$taxa}[0] = $add; }
        $taxa_add{$taxa}[1] = $add;
        $add = tell(BPOFILE);
      }
    }
    close(BPOFILE);
    $taxa_add{$taxa}[1] = $add;

    # point to global hash
    %taxa_bpo_index = %taxa_add;

    # export to disk
    if(!defined(nstore(\%taxa_bpo_index,$taxa_index_file)))
    {
      print "# failed while storing $taxa_index_file\n";
      unlink($taxa_index_file);
    }
  }

  print("# construct_taxa_indexes: number of taxa found = ".scalar(keys(%taxa_bpo_index))."\n");

#foreach $taxa (keys(%taxa_bpo_index)){ print "$taxa $taxa_bpo_index{$taxa}[0] to $taxa_bpo_index{$taxa}[1]\n"; }
}

# creates global hash index %blastquery
# considers only taxa in %only_this_taxa if requested, to save memory
# uses globals: %taxa_bpo_index,%gindex,@gindex2 (will populate them if necessary)
# Modified by Bruno Oct2015
sub construct_indexes
{
  my ($bpo_file,%only_this_taxa) = @_;

  my ($lastquery,$lasthit,$n_of_taxa,$q,$add) = ('',0,0);
  my ($doGindex,%taxon2name,@index,$taxon) = (0);

  $n_of_taxa = scalar(keys(%only_this_taxa));

  # sort taxa as they appear in $bpo_file, by checking global %taxa_bpo_index
  if($n_of_taxa)
  {
    $n_of_taxa = 0;
    foreach $taxon (sort {$taxa_bpo_index{$a}[0] <=> $taxa_bpo_index{$b}[0]} (keys(%taxa_bpo_index)))
    {
      next if(!$only_this_taxa{$taxon}); #print "$taxon $taxa_bpo_index{$taxon}[0]\n";
      $taxon2name{$n_of_taxa} = $taxon;
      $index[$n_of_taxa][0] = $taxa_bpo_index{$taxon}[0];
      $index[$n_of_taxa++][1] = $taxa_bpo_index{$taxon}[1];
    }

    %blastquery=(); # reset global variable
    if(!%gindex && !@gindex2){ $doGindex = 1 }
  }

  # parse BPO and leap as necessary to get to taxa of interest
  eval
  {
    $taxon=0;
    open (BPOFILE,$bpo_file) or die("# construct_indexes : cannot open $bpo_file file");
    flock(BPOFILE, 1);
    if($n_of_taxa){ seek(BPOFILE,$index[$taxon][0],0); }
    $add = tell(BPOFILE);
    while (<BPOFILE>)
    {
      #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      if(/^(\d+)/)
      {
        $q = $1;

        if($lastquery ne $q)
        {
          $blastquery{$q}[0] = $add;
			 
          if($lastquery)
          { 
            $blastquery{$lastquery}[1] = $lasthit; 
            $lasthit = 0; 
          }

          if($n_of_taxa && $add>$index[$taxon][1]) ## move to next taxon's block
          {
            $taxon++;
            last if($taxon==$n_of_taxa);
            seek(BPOFILE,$index[$taxon][0],0);
            $add = tell(BPOFILE); #print "mira: $add $index[$taxon][0] $index[$taxon-1][1]\n";
            next;
          }

          if($doGindex)
          {       
            # useful for cluster jobs, for which constructAllFasta was not invoked
            if(!defined($gindex{$taxon2name{$taxon}})){ $gindex{$taxon2name{$taxon}}[0] = $q }
            else
				    {
					    $gindex{$taxon2name{$taxon}}[1] = $q;
					    $gindex{$taxon2name{$taxon}}[2]++;
				    }
            $gindex2[$q] = $taxon2name{$taxon}; #print "$taxon $taxon2name{$taxon} $q\n";
          } 
        }

        $lasthit++;	
        $lastquery = $q;
        $add = tell(BPOFILE);
      }
    }
    close (BPOFILE);
    $blastquery{$lastquery}[1] = $lasthit;

    # point to global %blastquery
    $blastquery_ref = \%blastquery;
  };

  if($@){ die "# construct_indexes: failed: $@"; }

  printf("# number of file addresses/BLAST queries = %1.1e\n",scalar(keys(%blastquery)));

  open_bpofile($bpo_file); # leave open filehandle for future seek commands
} ## construct_indexes

# This subroutine is used to retrieve blast information from bpo file
# (blast parse out file), given the line_id.
# One Argument:
# 1. Reference Variable: reference to the address array
# 2. Number Variable: line_id of bpo file
# Modified by Bruno July2010
#sub getline_from_bpofile
#{
#  my ($lineid) = @_;

#  seek(BPOFILE,$address_ref->{$lineid},0); # used to be $address_ref->[$lineid-1]
#  my $line=<BPOFILE>;
#  my @bpout = split (";",$line);
#  chomp($bpout[8]); # $line=~s/[\r|\n]//g;
#  my ($pm,$pe) = ($bpout[5],0); #mantisa, exponente
#  if($bpout[5]=~/(\S+)e(\-\S+)/){ $pm=$1; $pe=$2; }

#  return ($bpout[0],$bpout[1],$bpout[2],$bpout[3],,$bpout[4],$pm,$pe,$bpout[6],$bpout[7],$bpout[8]); # this is why simspan goes from 7 to 8 !!!!
#} ## getline_from_bpofile                              ^^

# Three args:
# 1. file address
# 2. number of lines in block
# 3. number of columns to parse
# Edited on Oct2015 Bruno
# changed ; to \t
sub getblock_from_bpofile
{
  my ($add,$n_of_lines,$n_of_cols) = @_;

  my ($line,@block); 

  # find right line in file
  seek(BPOFILE,$add,0);

  if($n_of_cols)
  {	
	 while($n_of_lines > 0)
    {
	   $line=<BPOFILE>;
      chomp($line);
      push(@block,[split("\t",$line,$n_of_cols+1)]);
      $n_of_lines--;
    }
  }
  else
  {
	 while($n_of_lines > 0)
    {
      $line=<BPOFILE>;
      chomp($line);
      push(@block,[split("\t",$line)]);
      $n_of_lines--;
    }
  }

  return \@block;
}

# This subroutine is used to retrieve blast hits from bpo file
# given the query gene id and the subject gene id.
# If such information doesn't exist, zero is returned.
# Two Arguments:
# 1. String Variable: query gene id a
# 2. String Variable: subject gene id b
# Modified by Bruno Oct2015
sub blastqueryab
{
  my ($qa,$qb) = @_;

  if($blastquery_ref && defined($blastquery_ref->{$qa}))
  {
    my $ref_lines = getblock_from_bpofile($blastquery_ref->{$qa}[0],$blastquery_ref->{$qa}[1]);
    foreach my $line (@{$ref_lines})
    {
      if($line->[1] eq $qb){ return(@$line) }
    }
    return 0;
  }
  else{ return 0 }
} ## blastqueryab

# This subroutine is used to generate an input file (matrix file)
# for MCL and index file.
# uses global variables $matrix_file, @graph, %weight, $last_graph_item
# Modified by Bruno May2012
sub write_matrix_index
{
  my ($i,$p,$m) = (1);
  my $size = $last_graph_item+1; #scalar(@graph) fails with BerkeleyDB tie

  open (MTX,">$matrix_file") or die("# write_matrix_index: cannot write to file $matrix_file");
  print MTX "(mclheader\nmcltype matrix\ndimensions ".$size."x".$size."\n)\n\n(mclmatrix\nbegin\n\n";
  foreach $p (1 .. $last_graph_item)
  {
    print MTX "$p     ";
    my @nodes = sort split(/,/,$graph[$p]);
    foreach $m (@nodes){ print MTX $m.':'.$weight{$p.' '.$m}.' ' }
    print MTX "\$\n";
  }
  print MTX ")\n\n";
  close (MTX); #print("\nMatrix($size X $size) file $matrix_file generated\n");
}

# inspired by original mcl_backindex
#uses globals: @gindex2,
# created by Bruno, last changed Oct2015
sub parse_MCL_clusters
{
  my ($mcl_file) = @_;

  my (%orthologues,%orth_taxa);
  my ($last,$lastline,$cluster,$orth,$ref,$taxon,@mcl) = (0,'');

  open (MCL,$mcl_file) or die("# parse_MCL_clusters: can't open $mcl_file");
  while (<MCL>)
  {
    chomp; s/\$$//; last if (/^\)/ && scalar(@mcl));
    if (/^(\d+)\s+(.*)\$/ || /^(\d+)\s+(.*)/)
    {
      $mcl[$last]=$lastline;
      $last=$1;$lastline=$2;
    }
    elsif (/^\s+/) {$lastline.=$_;}
  }
  $mcl[$last]=$lastline;
  close (MCL);

  foreach $cluster (@mcl)
  {
    my (@cluster) = (split(/\s+/,$cluster));

    # record taxa composition as well
    $ref = $cluster[0];
    foreach $orth (@cluster)
    {
      $taxon = $gindex2[$orth];
      $orth_taxa{$ref}{$taxon}++;
    }

    # record orthologus cluster
    shift(@cluster);
    foreach $orth (@cluster){ push(@{$orthologues{$ref}},$orth); }
  }

  return (\%orthologues,\%orth_taxa);
} ## parse_MCL_clusters

sub get_makeIsoform_outfilename
{
  my ($taxon) = @_;
  return $TMP_DIR.'isoforms_'.$taxon;
}

# Modified from makeInparalog, uses global variables $bpo_file,@gindex2,$TMP_DIR
# assumes sequences are pre-sorted from largest to smallest
# Bruno Oct2015, April2016
sub makeIsoform
{
  my ($saveRAM,$taxon,$evalue_cutoff,$min_overlap,$same_best_hit,$force_parsing) = @_;

  my $refh_redundant_isoforms;
  my $iso_file = get_makeIsoform_outfilename($taxon);

  # 0) re-use previously calculated isoforms if possible
  if(-s $iso_file && !$force_parsing)
  {

    my ($ref_data_structure);

    eval { $ref_data_structure = retrieve($iso_file) };
    if(!$@) #  compatible previous binary $inpara_file (avoid trans-arch errors)
    {
      if(defined($ref_data_structure))
      {
        $refh_redundant_isoforms = $ref_data_structure;

        printf("# %s : %d sequences\n",$taxon,scalar(keys(%$refh_redundant_isoforms)));

        return ($refh_redundant_isoforms);
      }
    }
  }

  # 1) otherwise find isoforms by parsing blast output (default behaviour)
  my ($qid,$sid,$ev,$pi);
  my ($overlap,$concat_residues,$no_overhang,$line);
  my ($qlength,$slength,$span,$qstart,$qend,$sstart,$send);

  if($saveRAM){ construct_indexes($bpo_file,($taxon=>1)) }

  # 1.1) find overlapping isoforms
  foreach $qid ($gindex{$taxon}[0] .. $gindex{$taxon}[1])
  {
    next if($refh_redundant_isoforms->{$qid} || !defined($blastquery{$qid}));

    # assumes in-species hits are first
    my $ref_lines = getblock_from_bpofile($blastquery{$qid}[0],$blastquery{$qid}[1]);

    foreach $line (@$ref_lines)
    {
      #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits

		#1	5406	0.0	99	13	100	15307	2119	1:13190-15307:1-2119.	3906
		#1	54494	0.0	99	99	99	15307	15270	1:1-15265:15269-5.	2.813e+04
      ($sid,$ev,$pi,$qlength,$slength,$span) = @$line[1,2,3,6,7,8];

      next if($sid <= $qid || $refh_redundant_isoforms->{$sid});
      next if($pi < 100); # allow only perfect matches
      next if($span =~ /\.2:\d+/); # allow only matches with a single,contiguous hsp
      last if($ev > $evalue_cutoff);
      next if(!defined($gindex2[$sid]) || $gindex2[$sid] ne $taxon); # skip hits from other species

      # check hits that match 5' or '3 ends with no overhangs
      $no_overhang = 0;
      if($span =~ m/^1:(\d+)\-(\d+)\:(\d+)\-(\d+)/)
      {
        ($qstart,$qend,$sstart,$send) = ($1,$2,$3,$4);

        # q -------     or  --------
        # s      -----        ------
        if($qend == $qlength && $sstart == 1)
        {
          $no_overhang = 1;
          $concat_residues = $send;
        }

        # q   -------  or    -----
        # s ----           -------
        elsif($qstart == 1 && $send == $slength)
        {
          $no_overhang = 1;
          $concat_residues = $qend;
        }
      }

      if($no_overhang)
      {
        next if($concat_residues < $min_overlap);
        #print "noover $qid $sid $concat_residues,$pi $span $qlength $slength\n";
      }
      else
      {
        ($overlap,$concat_residues) = simspan_hsps_overlap($qlength,$slength,$span);
        next if($overlap < 99 || $concat_residues < $min_overlap);
        #print "overla $qid $sid $overlap,$concat_residues,$pi $span $qlength $slength\n";
      }

      $refh_redundant_isoforms->{$sid} = $qid;
    }
  }#printf("# %d overlapping isoforms\n",scalar(keys(%$refh_redundant_isoforms)));

  # 1.2) optionally merge transcript fragments, with insufficient or null overlap, matching the same best hit
  # as before, conserve the longest one (not tested 12012015)
  if($same_best_hit)
  {
    my ($n_of_besthit_isoforms,$refid,%besthit) = (0);
    foreach $qid ($gindex{$taxon}[0] .. $gindex{$taxon}[1])
    {
      next if($refh_redundant_isoforms->{$qid} || !defined($blastquery{$qid}));

      my $ref_lines = getblock_from_bpofile($blastquery{$qid}[0],$blastquery{$qid}[1],3);
      foreach $line (@$ref_lines)
      {
        #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
        ($sid,$ev) = @$line[1,2];

        next if($refh_redundant_isoforms->{$sid} || $sid eq $qid || $gindex2[$sid] eq $taxon);
        last if($ev > $evalue_cutoff);
        push(@{$besthit{$sid}},$qid); # keep adding qids, which are by definition length-sorted
        last; # take only best hit
      }
    }

    # take longest transcripts matching same best hit
    foreach $sid (keys(%besthit))
    {
      next if(scalar(@{$besthit{$sid}})<2);
      my @hits = @{$besthit{$sid}};
      $refid = shift(@hits);
      foreach $sid (1 .. $#hits)
      {
        $refh_redundant_isoforms->{$sid} = $refid;
        $n_of_besthit_isoforms++;
      }
    }

    printf("# %d isoforms sharing best hit\n",$n_of_besthit_isoforms);
  }

  printf("# %s : %d sequences\n",$taxon,scalar(keys(%$refh_redundant_isoforms)));

  ## save hash to file with Storable::store function
  if(!defined(nstore($refh_redundant_isoforms,$iso_file)))
  {
    print "# makeIsoform: failed while storing isoforms to $iso_file\n";
    unlink($iso_file);
  }

  return $refh_redundant_isoforms;
}

sub get_makeInparalog_outfilename
{
  my ($taxon) = @_;
  return $TMP_DIR.'inparalogues_'.$taxon;
}

# This subroutine is an important part of OrthoMCL, used to look for inparalogs, defined as
# reciprocal better hits encoded in the same genome. It is used also to look for BDBHs.
# 8 arguments:
# 1. String: flag used to request indexes and save RAM
# 2. String Variable: Taxon name
# 3,4,5,6. String: E-value, Identity, Coverage and Neighbor correlation cut-offs for potential homologous
# 7 flag that says whether inparalogous blast hits are already sorted, otherwise they need to be sorted
# 8 flag to force $bpo_file parsing, which can otherwise be avoided if file exists
# 9 flag to indicate that shortest sequence should be used to calculate coverage
# Modified by Bruno Oct2011,Oct2015 uses global variables $bpo_file,%gindex,@gindex2,$TMP_DIR
sub makeInparalog
{
  my ($saveRAM,$taxon,$evalue_cutoff,$pi_cutoff,
    $pmatch_cutoff,$neighbor_corr_cutoff,
    $inparalogues_are_sorted,$force_parsing,$use_short_sequence) = @_;

  if(!$use_short_sequence){ $use_short_sequence = 0 }

  my $inpara_file = get_makeInparalog_outfilename($taxon);
  my ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);

  if(-s $inpara_file && !$force_parsing)
  {

    # try to re-use previously calculated inparalogues
    my ($ref_data_structure);

    eval { $ref_data_structure = retrieve($inpara_file) };
    if(!$@) #  compatible previous binary $inpara_file (avoid trans-arch errors)
    {
      if(defined($ref_data_structure))
      {
        $ref_hash_edges = $ref_data_structure->{'edges'};
        $ref_hash_weights = $ref_data_structure->{'weights'};
        $avgw = $ref_data_structure->{'average_weight'};
        $sumw = $ref_data_structure->{'sum_weights'};
        $count = $ref_data_structure->{'count'};

        printf("# %d sequences\n",scalar(keys %$ref_hash_weights)/2);

        return ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);
      }
    }
  }

  # default behaviour, calculate best hits by parsing blast output
  my (%inbest,%pvalue,$qid,$sid,$ev,$pi,$line);
  my ($qcov,$scov,%data_structure);

  if($saveRAM){ construct_indexes($bpo_file,($taxon=>1)) }

  foreach $qid ($gindex{$taxon}[0] .. $gindex{$taxon}[1])
  {
    next if(!defined($blastquery{$qid}));
    my $ref_lines;

    if($inparalogues_are_sorted){ $ref_lines = getblock_from_bpofile($blastquery{$qid}[0],$blastquery{$qid}[1],6); }
    else
    {
      # resort $qid hits, putting them ahead of equal E-value hits from other species
      $ref_lines = sort_taxon_hits($blastquery{$qid}[0],$blastquery{$qid}[1],$taxon);
    }

    foreach $line (@$ref_lines)
    {
		  #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      ($sid,$ev,$pi,$qcov,$scov) = @$line[1,2,3,4,5];

      next if($sid eq $qid); # ignore autohits
      last if(!defined($gindex2[$sid]) || $gindex2[$sid] ne $taxon); # best hit from other species
      last if($ev > $evalue_cutoff || $pi < $pi_cutoff);
		  next if($pmatch_cutoff && ( 
			($use_short_sequence && $qcov < $pmatch_cutoff && $scov < $pmatch_cutoff) ||
			(!$use_short_sequence && ($qcov < $pmatch_cutoff || $scov < $pmatch_cutoff)) ));
      
      push(@{$inbest{$qid}}, $sid);
      $pvalue{$qid.' '.$sid} = $ev;  # pareja de secuencias (query,subject) con su e-value
    }
  }

  #print("# ".scalar(keys(%inbest))." sequences have best hits within species\n");

  ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count) = matrix(\%inbest, \%pvalue, $neighbor_corr_cutoff);

  # group these variables in a single hash for Storable::store function
  %data_structure = ( 'edges'=>$ref_hash_edges,'weights'=>$ref_hash_weights,
    'average_weight'=>$avgw,'sum_weights'=>$sumw,'count'=>$count );

  if(!defined(nstore(\%data_structure,$inpara_file)))
  {
    print "# makeInparalog: failed while storing inparalogues to $inpara_file\n";
    unlink($inpara_file);
  }

  return ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);
}


## This subroutine is a variation of makeInparalog for clustered sequences, which can then be used 
## by OMCL and BDBH algorithms. It calculates 
## Arguments:
## 1. String: flag used to request indexes and save RAM
## 2. String Variable: Taxon name
## 3. Hash reference with sequence ids (values) of clusters (keys)
## 4. flag to force cluster parsing, which can otherwise be avoided if file exists
## uses global variable $bpo_file,%gindex,@gindex2,$TMP_DIR
## mar2015 Bruno  
#sub makeClusterInparalog
#{
#  my ($saveRAM,$taxon,$ref_clusters,$force_parsing) = @_;
#
#  my $inpara_file = get_makeInparalog_outfilename($taxon);
#  my ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);
#
#  if(-s $inpara_file && !$force_parsing)
#  {
#    # try to re-use previously calculated inparalogues
#    my ($ref_data_structure);
#
#    eval { $ref_data_structure = retrieve($inpara_file) };
#    if(!$@) #  compatible previous binary $inpara_file (avoid trans-arch errors)
#    {
#      if(defined($ref_data_structure))
#      {
#        $ref_hash_edges = $ref_data_structure->{'edges'};
#        $ref_hash_weights = $ref_data_structure->{'weights'};
#        $avgw = $ref_data_structure->{'average_weight'};
#        $sumw = $ref_data_structure->{'sum_weights'};
#        $count = $ref_data_structure->{'count'};
#
#        printf("# %d sequences\n",scalar(keys %$ref_hash_weights)/2);
#
#        return ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);
#      }
#    }
#  }
#
#  # explicitly call inparalogues as members of the same cluster
#  my ($cl,$qid,$sid,$pm,$pe,$w,$pv1,$pv2,$sStart,$sEnd,$line);
#  my (%data_structure,%edges,%weights,%self_evalues);
#  
#  
#  if($saveRAM){ construct_indexes($bpo_file,($taxon=>1)) }
#
#  foreach $qid (@{$gindex{$taxon}})
#  {
#    next if(!defined($blastquery{$qid}));
#    ($sStart,$sEnd) = split (";",$blastquery{$qid});
#    my $ref_lines = getblock_from_bpofile($sStart,$sEnd-$sStart+1); 
#    foreach $line (@$ref_lines)
#    {
#      #1;1;521;1;521;0.0;100;1:1-521:1-521.;1025
#      #2;1;521;1058;522;2e-173;60;1:14-515:14-515.;610
#      ($sid,$pm,$pe) = @$line[3,5,6];
#      if($sid eq $qid)
#      {
#        $self_evalues{$qid} = $pm.'e'.$pe;
#        last;
#      }
#    }
#  }
#  
#  $count = 0;
#  foreach $cl (sort {$a<=>$b} keys(%$ref_clusters))
#  { 
#    my @ids = @{$ref_clusters->{$cl}}; #1991,1992,1993,1994 print "$cl ".join(',',@ids)."\n";
#    
#    foreach $qid (0 .. $#ids-1)
#    {
#      foreach $sid ($qid+1 .. $#ids)
#      {
#        #1991: 1992,1993,1994
#        #1992: 1993,1994
#        #1993: 1994
#        
#        push( @{$edges{$ids[$qid]}}, $ids[$sid] );
#        push( @{$edges{$ids[$sid]}}, $ids[$qid] );
#        
#        $pv1 = $self_evalues{$ids[$qid]};
#        $pv2 = $self_evalues{$ids[$sid]};
#        
#        if($pv1 == 0) { $pv1 = $MAX_WEIGHT_DEFAULT }
#        else{ $pv1 = -log($pv1)/log(10) }
#        if($pv2 == 0) { $pv2 = $MAX_WEIGHT_DEFAULT }
#        else{ $pv2 = -log($pv1)/log(10) }
#        
#        $w = ($pv1+$pv2)/2;
#        $sumw += $w;
#        $count++;
#        
#        $weights{$ids[$qid].' '.$ids[$sid]} = sprintf("%.3f", $pv1);
#        $weights{$ids[$sid].' '.$ids[$qid]} = sprintf("%.3f", $pv2);
#      }
#    }
#  }
#  
#  $avgw = 'N/A';
#  if($count){ $avgw = $sumw/$count }
#  
#  print("# $count sequences\n"); 
#
#  # group these variables in a single hash for Storable::store function
#  %data_structure = ( 'edges'=>\%edges,'weights'=>\%weights,
#    'average_weight'=>$avgw,'sum_weights'=>$sumw,'count'=>$count );
#
#  if(!defined(nstore(\%data_structure,$inpara_file)))
#  {
#    print "# makeClusterInparalog: failed while storing inparalogues to $inpara_file\n";
#    unlink($inpara_file);
#  }
#
#  return (\%edges,\%weights,$avgw,$sumw,$count);
#}


sub get_makeOrtholog_outfilename
{
  my ($options,$taxon1,$taxon2) = @_;
  return $TMP_DIR.'orthologues_'.$options.'_'.$taxon1.'_'.$taxon2;
}

# Up to 13 arguments
# 1. String: flag used to request indexes and save RAM
# 2. String Variable: Taxon name A, assumed to be reference
# 3. String Variable: Taxon name B
# 4. String Variable: flag stating that only best hit should be taken
# 5,6,7,8. String: E-value, Identity, Coverage and Neighbor correlation cut-offs for potential homologous
# 9: (Optional) String: flag meaning that no reference genome should be considered when calling matrix()
# 10: (Optional) String: flag for requesting that orthologs should have identical Pfam domain strings, it might be a path combined with 1
# 11: (Optional) String: flag for requesting explicitly that blast output needs to be reparsed
# 12,13: (Optional) Hash reference: hash of inparalogue seeds, filled in
# 14: flag to indicate that shortest sequence should be used to calculate coverage
# Modified by Bruno Nov2011,Jan2015
sub makeOrtholog
{
  my ($saveRAM,$ta,$tb,$only_best_match,$evalue_cutoff,$pi_cutoff,$pmatch_cutoff,
    $min_neighbor_cutoff,$noreference,$checkPFAMdoms,$force_parsing,
    $rhash_inparala,$rhash_inparalb,$use_short_sequence) = @_;

  my $pfam_file_path = '';
  if($checkPFAMdoms)
  {
    $pfam_file_path = $checkPFAMdoms; # in case a path is passed for cluster jobs
    $checkPFAMdoms = 1;
  }

  if(!$use_short_sequence){ $use_short_sequence = 0 }
  else{ $use_short_sequence = 1 }

  my $orth_file = get_makeOrtholog_outfilename($only_best_match.$noreference.$checkPFAMdoms,$ta,$tb);

  my ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);

  if(-s $orth_file && !$force_parsing)
  {

    # try to re-use previously calculated orthologues
    my ($ref_data_structure);

    eval { $ref_data_structure = retrieve($orth_file) };
    if(!$@) #  compatible previous binary $orth_file (avoid trans-arch errors)
    {
      if(defined($ref_data_structure))
      {
        $ref_hash_edges = $ref_data_structure->{'edges'};
        $ref_hash_weights = $ref_data_structure->{'weights'};
        $avgw = $ref_data_structure->{'average_weight'};
        $sumw = $ref_data_structure->{'sum_weights'};
        $count = $ref_data_structure->{'count'};

        printf("# %d sequences\n",scalar(keys %$ref_hash_weights)/2);

        return ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);
      }
    }
  }

  # default behaviour, calculate best hits by parsing blast output
  my (%best,%sim,%pvalue,$pfamq,$pfams,$sid,$ev,$pi,%data_structure);
  my ($lastev,$hit_id,$qid,$q,$line,$qcov,$scov);

  if($saveRAM)
  {
    construct_indexes($bpo_file,($ta=>1,$tb=>1));
    if($checkPFAMdoms){ construct_Pfam_hash($pfam_file_path) }
  }

  foreach $qid ($gindex{$ta}[0] .. $gindex{$ta}[1])
  {
    next if(!defined($blastquery{$qid}));
    ($lastev,$hit_id) = ('',0);

    my $ref_lines = getblock_from_bpofile($blastquery{$qid}[0],$blastquery{$qid}[1],6);
    foreach $line (@{$ref_lines})
    {
      #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      ($sid,$ev,$pi,$qcov,$scov) = @$line[1,2,3,4,5];

      if(defined($gindex2[$sid]))
      {
        if($gindex2[$sid] eq $tb)
        {
          $hit_id++;

          # try to collapse inparologous matches
          $qid = $rhash_inparala->{$qid} || $qid;
          $sid = $rhash_inparalb->{$sid} || $sid;

          if($hit_id==1)
          {
            push(@{$sim{$qid}},[$sid,$ev,$pi,$qcov,$scov]);
            $lastev=$ev;
            #print "$qid $sid,$pm,$pe,$pi,$qlength,$slength,$span\n";
          }
          else
          {
            if($lastev eq $ev)
            {
              push(@{$sim{$qid}},[$sid,$ev,$pi,$qcov,$scov]);
            }
            else{last}
          }
        }
      } #else {print("$sid gindex2 not defined; lineid: $_\n");}
    }
  }
  foreach $qid ($gindex{$tb}[0] .. $gindex{$tb}[1])
  {
    next if(!defined($blastquery{$qid}));
    ($lastev,$hit_id) = ('',0);

    my $ref_lines = getblock_from_bpofile($blastquery{$qid}[0],$blastquery{$qid}[1],6);
    foreach $line (@{$ref_lines})
    {
		#$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      ($sid,$ev,$pi,$qcov,$scov) = @$line[1,2,3,4,5];
      if (defined $gindex2[$sid])
      {
        if ($gindex2[$sid] eq $ta)
        {
          $hit_id++;

          # try to collapse inparologous matches
          $qid = $rhash_inparalb->{$qid} || $qid;
          $sid = $rhash_inparala->{$sid} || $sid;

          if ($hit_id==1)
          {
            push(@{$sim{$qid}},[$sid,$ev,$pi,$qcov,$scov]);
            $lastev=$ev;
            #print "$qid $sid,$pm,$pe,$pi,$qlength,$slength,$span\n";
          }
          else
          {
            if($lastev eq $ev)
            {
              push(@{$sim{$qid}},[$sid,$ev,$pi,$qcov,$scov]);
            }
            else{last}
          }
        }
      } #else {print("$sid gindex2 not defined; lineid: $_\n");}
    }
  }

  foreach $q (keys(%sim))
  {
    foreach (@{$sim{$q}})
    {
      ($sid,$ev,$pi,$qcov,$scov) = @$_;
      if($checkPFAMdoms)
      {
        $pfamq = $pfam_hash{$q} || '';
        $pfams = $pfam_hash{$sid} || '';
        if($pfamq ne $pfams)
        {
          # currently requires identical Pfam annotation, however other options could be considered
          #print "# skip BDBH with different Pfam annotations: $q|$pfamq $bla[0]|$pfams\n";
          next;
        }
      }

      # skip matches of low sequence identity
      next if($pi < $pi_cutoff);

      # skip matches of poor two-ways coverage
      # check coverage on the other BLAST direction as alignments are frequently
      # different and might well improve the coverage (example lspA among Gammaproteobacteria)

		  if($pmatch_cutoff && (
         ($use_short_sequence && $qcov < $pmatch_cutoff && $scov < $pmatch_cutoff) ||
         (!$use_short_sequence && ($qcov < $pmatch_cutoff || $scov < $pmatch_cutoff)) ))
      {
        next if($noreference);
        my $pmatch2way = 0;
        foreach my $hit2way (@{$sim{$sid}})
        {
          next if($hit2way->[0] ne $q);
	
			 if( ($use_short_sequence && ($hit2way->[3] >= $pmatch_cutoff || $hit2way->[4] >= $pmatch_cutoff)) ||
				(!$use_short_sequence && ($hit2way->[3] >= $pmatch_cutoff && $hit2way->[4] >= $pmatch_cutoff)) )
			 {
				$pmatch2way = 1 
			 }
          last;
        }
        next if(!$pmatch2way);
      }

      # skips matches of poor two-ways E-value
      # short sequences might produce poor Evalues (example rpmH among Gammaproteobacteria,
      # 46aa with 3e-05, just over default cutoff)
      if($ev > $evalue_cutoff)
      {
        next if($noreference);
        my $Evalue2way = 0;
        foreach my $hit2way (@{$sim{$sid}})
        {
          next if($hit2way->[0] ne $q);
          if($hit2way->[1] <= $evalue_cutoff){ $Evalue2way = 1 }
          last;
        }
        next if(!$Evalue2way);
      }

      push(@{$best{$q}}, $sid);
      $pvalue{$q.' '.$sid} = $ev; #print "$q $sid\n" if($q == 1119 || $sid == 1119); # 1119 , 77 debug

      last if($only_best_match); # you might find >1 bdbhs otherwise, provided they pass the 3 cutoffs
    }
  }

  #print("# ".scalar(keys(%best))." sequences have best hits from other species\n");
  if($noreference){ ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count) = matrix(\%best, \%pvalue, $min_neighbor_cutoff); }
  else{ ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count) = matrix(\%best, \%pvalue, $min_neighbor_cutoff, $ta); }

  # group these variables in a single hash for Storable::store function
  %data_structure = ( 'edges'=>$ref_hash_edges,'weights'=>$ref_hash_weights,
    'average_weight'=>$avgw,'sum_weights'=>$sumw,'count'=>$count );

  if(!defined(nstore(\%data_structure,$orth_file)))
  {
    print "# makeOrtholog: failed while storing orthologues to $orth_file\n";
    unlink($orth_file);
  }

  return ($ref_hash_edges,$ref_hash_weights,$avgw,$sumw,$count);

} ## makeOrtholog

sub get_makeHomolog_outfilename
{
  my ($taxon1,$taxon2) = @_;
  return $TMP_DIR.'homologues_'.$taxon1.'_'.$taxon2;
}

# looks for best homolog in B for sequences in A, used to calculate pangenome
# 10 arguments
# 1. String: flag used to request indexes and save RAM
# 2. String Variable: Taxon name A, assumed to be reference
# 3. String Variable: Taxon name B
# 4. String Variable: flag stating that only best hit should be taken
# 5,6,7,8. String: E-value, Identity, Coverage and Neighbor correlation cut-offs for potential homologous
# 9 flag: force blast out parsing
# 10: flag to indicate that coverage should be calculated with respect to shortest sequence
# 11 flag to indicate that homologues should be searched in both directions
# Written by Bruno in sept2010,Jan2015
sub makeHomolog
{
  my ($saveRAM,$ta,$tb,$evalue_cutoff,$pi_cutoff,$pmatch_cutoff,$force_parsing,$use_short_sequence,$bidi) = @_;

  if(!$use_short_sequence){ $use_short_sequence = 0 }

  my $homol_file = get_makeHomolog_outfilename($ta,$tb);
  if(-s $homol_file && !$force_parsing)
  {
    # try to re-use previously calculated homologues
    my ($ref_data_structure);

    eval{ $ref_data_structure = retrieve($homol_file) };
    if(!$@)
    {
      if(defined($ref_data_structure))
      {
        print("# ".scalar(keys(%$ref_data_structure))." sequences \n");#$homol_file)\n");
        return $ref_data_structure;
      }
    }
  }

  # default behaviour, calculate best hits by parsing blast output
  my (%homologues,$qid,$sid,$ev,$pi,$line,$qcov,$scov);

  if($saveRAM){ construct_indexes($bpo_file,($ta=>1,$tb=>1)) }

  foreach $qid ($gindex{$ta}[0] .. $gindex{$ta}[1]) # search for homologs in B for A queries
  {
    next if(!defined($blastquery{$qid}));

    my $ref_lines = getblock_from_bpofile($blastquery{$qid}[0],$blastquery{$qid}[1],6);
    foreach $line (@{$ref_lines})
    {
      #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      ($sid,$ev,$pi,$qcov,$scov) = @$line[1,2,3,4,5];
      next if(!defined($gindex2[$sid]) || $gindex2[$sid] ne $tb);
      last if($ev > $evalue_cutoff);
      next if($pi < $pi_cutoff);

		  next if($pmatch_cutoff && (
         ($use_short_sequence && $qcov < $pmatch_cutoff && $scov < $pmatch_cutoff) ||
         (!$use_short_sequence && ($qcov < $pmatch_cutoff || $scov < $pmatch_cutoff)) ));
      $homologues{$qid} = $sid;
      last;
    }
  }

  if($bidi) # search for homologs in A for B queries in case there are unidirectional hits
  {
#foreach $qid (@{$gindex{$tb}}){
#	next if(!defined($blastquery{$qid}));
#        ($sStart,$sEnd) = split(";",$blastquery{$qid});
#       my $ref_lines = getblock_from_bpofile($sStart,$sEnd-$sStart+1);
#       foreach $line (@{$ref_lines}){
#		($qlength,$sid,$slength,$pm,$pe,$pi,$span) = @$line[2,3,4,5,6,7,8];
#		last if($homologues{$qid});
#		if(defined $gindex2{$sid}) {
#			next if($gindex2{$sid} ne $ta);
#			last if($pm.'e'.$pe > $evalue_cutoff);
#                        next if($pi < $pi_cutoff); # es el corte más restrictivo
#                        next if($pmatch_cutoff && simspan_hsps($qlength,$slength,$span) < $pmatch_cutoff);
#			$homologues{$sid} = $qid;
#			last;
#		}
#	}
#}
  }

  print("# ".scalar(keys(%homologues))." sequences\n"); # have homologues from the other species\n");

  # store hash in file with Storable::store
  if(!defined(nstore(\%homologues,$homol_file)))
  {
    print "# failed while storing homologues to $homol_file\n";
    unlink($homol_file);
  }

  return \%homologues ;
} ## makeHomolog

# Defines clusters of lineage specific expansions (inparalogues).
# Takes 1 arg, a reference to a %edge hash returned by matrix()
# Returns a reference to a hash storing the clusters
# Created by Bruno Jun2010
sub cluster_lineage_expansions
{
  my ($ref_hash_inparalogues) = @_;

  my (%inparalog_seed);
  foreach my $gene (keys(%{$ref_hash_inparalogues}))
  {
    # take as seed or representative inparalog with smallest id
    my $seed = (sort {$a<=>$b} (@{$ref_hash_inparalogues->{$gene}}))[0];
    next if($seed > $gene); #print "# $gene $rep\n";
    $inparalog_seed{$gene} = $seed;
  }

  return (\%inparalog_seed);
}

# This subroutine is used to choose two-way hits among one-way hits (best
# hits between two species or better hits within one species),
# calculate the weight between two nodes (minus logarithm of the p-value,
# or $MAX_WEIGHT_DEFAULT for p-value 0 ), and calculate average
# weight among all inparalogs within one species or all orthologs between
# two species. (Weighting process takes place in the main script)
# => %edge contiene parejas llave,valor de bdbhs
# Two Arguments:
# 1. Reference Variable: reference to a hash which stores all the possible
#    gene pairs (one-way best hit, or better hit).
# 2. Reference Variable: reference to a hash which stores the pvalue for
#    the gene pairs.
# uses globals: @gindex2
# Modified by Bruno in Aug2008,Oct2015
sub matrix
{
  my ($ref_best,$ref_pvalue,$min_neighbor_corr,$reference_taxon) = @_;

  my ($count,$sumw,$flag,$pv1,$pv2,$w,%edge, %weight) = (0,0,0);

  foreach my $query (sort(keys(%$ref_best)))
  {
    foreach my $subject (@{$ref_best->{$query}})
    {
      next if($weight{$query.' '.$subject});
      $flag = 0;

      # check this is really a bidirectional best hit
      foreach my $q (@{$ref_best->{$subject}}) { if($q eq $query) { $flag = 1; } }

      # check these sequences have enough neighbor correlation (NC)
      if($flag && $min_neighbor_corr &&
        neighbor_correlation($query,$subject) < $min_neighbor_corr){ $flag = 0 }

      if($flag == 1)
      {
        if($reference_taxon)
        {
          if($gindex2[$query] eq $reference_taxon){ push (@{$edge{$query}}, $subject); }
          else{ push (@{$edge{$subject}}, $query); }
        }
        else # default behaviour
        {
          push (@{$edge{$query}}, $subject);
          push (@{$edge{$subject}}, $query);
        }

        #use -logP as weights and treat P=0 as -logP=$maximum_weight (DEFAULT=300)
        if($ref_pvalue->{$query.' '.$subject} == 0) { $pv1 = $MAX_WEIGHT_DEFAULT; }
        else { $pv1 = -log($ref_pvalue->{$query.' '.$subject})/log(10); }
        if($ref_pvalue->{$subject.' '.$query} == 0) { $pv2 = $MAX_WEIGHT_DEFAULT; }
        else { $pv2 = -log($ref_pvalue->{$subject.' '.$query})/log(10); }

        $w = ($pv1+$pv2)/2;
        $sumw += $w;
        $count++;

        # use averaged score as edge weight
        $weight{$query.' '.$subject} = sprintf("%.3f", $w);
        $weight{$subject.' '.$query} = sprintf("%.3f", $w);
      }
    }
  }
  my $avgw = 'N/A';
  if ($count) { $avgw = $sumw/$count; }
  my $no_tmp = scalar(keys(%weight))/2;
  print("# $no_tmp sequences\n"); # pairs were identified as Reciprocal Best Hits\n");

  return (\%edge, \%weight, $avgw, $sumw, $count);
} ## matrix

# This subroutine is used by the subroutine makeInparalog, rearranging
# Evalue-tied blast hits so that the hits from a specific taxon are moved
# higher than hits from other species with same E-value
# Three Arguments:
# 1. Number Variable: address in bpo file 
# 2. Number Variable: last hit counting from address
# 3. String Variable: taxon
# Modified by Bruno Oct2015
sub sort_taxon_hits
{
  my ($add,$nhits,$taxon)=@_;

  my (@sorted_simid,@tmp,$sid,$ev,$lastev,$line);
  my $ref_lines = getblock_from_bpofile($add,$nhits,3);
  $line = $ref_lines->[0];
  $lastev = @$line[2];

  foreach $line (@{$ref_lines})
  {
    ($sid,$ev) = @$line[1,2];
    if($lastev eq $ev)
    {
      if ($gindex2[$sid] eq $taxon){	push (@sorted_simid,$line); }
      else { push (@tmp,$line); }
    }
    else
    {
      if(scalar(@tmp)>0){ push (@sorted_simid,@tmp); @tmp=(); }
      if($gindex2[$sid] eq $taxon){ push (@sorted_simid,$line); }
      else{ push (@tmp,$line); }
    }
    $lastev=$ev;
  }
  if (scalar(@tmp)>0) {push (@sorted_simid,@tmp);}

  return \@sorted_simid;
}

sub alignment_coords
{
  my ($hsps) = @_;

  my ($qstart,$qend,$sstart,$send) = (10000000000,0,10000000000,0);
  my ($q1,$q2,$s1,$s2);

  foreach (split ('\.',$hsps)) #$simspan .= "$n_of_hsp:$pQstart-$pQend:$pSstart-$pSend.";
  {
    if (/\d+\:(\d+)\-(\d+)\:(\d+)\-(\d+)/)
    {
      ($q1,$q2,$s1,$s2) = ($1,$2,$3,$4);

      if($q1 < $qstart){ $qstart = $q1 }
      if($s1 < $sstart){ $sstart = $s1 }
      if($q2 > $qend){ $qend = $q2 }
      if($s2 > $send){ $send = $s2 }
    }
  }

  return ($qstart,$qend,$sstart,$send);
} ## alignment_coords

# 3 args: $lengthq,$lengths,$hsps
# returns 1) % overlap with respect to short sequence and 2) number of overlapping residues
# Bruno Jan2015
sub simspan_hsps_overlap
{
  my ($lengthq,$lengths,$hsps) = @_;

  my (%sub_start, %sub_length, %query_start, %query_length, $overlap);
  my ($sub_segment_first,$sub_segment_last,$query_segment_first,$query_segment_last) = ('','','','');

  my @hsp=split ('\.',$hsps);#"1:10-22:2-14.2:25-75:15-65";
  foreach (@hsp)
  {
    if (/(\d+)\:(\d+)\-(\d+)\:(\d+)\-(\d+)/)
    {
      if($5>$4)
      {
        $sub_start{$1}=$4;
        $sub_length{$1}=$5-$4+1;

        if(!$sub_segment_first || $4 < $sub_segment_first){ $sub_segment_first = $4; }
        if(!$sub_segment_last  || $5 > $sub_segment_last){ $sub_segment_last = $5; }
      }
      else # subject in reverse strand
      {
        $sub_start{$1}=$5;
        $sub_length{$1}=$4-$5+1;

        if(!$sub_segment_first || $5 < $sub_segment_first){ $sub_segment_first = $5; }
        if(!$sub_segment_last  || $4 > $sub_segment_last){ $sub_segment_last = $4; }
      }

      $query_start{$1}=$2;
      $query_length{$1}=$3-$2+1;

      if(!$query_segment_first || $2 < $query_segment_first){ $query_segment_first = $2; }
      if(!$query_segment_last  || $3 > $query_segment_last){ $query_segment_last = $3; }
    }
  }

  if($lengthq > $lengths)
  {
    $overlap = matchlen(\%sub_start,\%sub_length);
    return ( 100*$overlap/$lengths , $overlap );
  }
  else
  {
    $overlap = matchlen(\%query_start,\%query_length);
    return ( 100*$overlap/$lengthq , $overlap );
  }
}

# 3 args: $lengthq,$lengths,$hsps
# Added by Bruno Aug2008, modified Jan2015
sub simspan_hsps
{
  my ($lengthq,$lengths,$hsps,$use_short_sequence) = @_;

  my (%sub_start, %sub_length, %query_start, %query_length);
  my @hsp=split ('\.',$hsps); #my $max_segment;
  my ($sub_segment_first,$sub_segment_last,$query_segment_first,$query_segment_last) = ('','','','');

  #$simspan .= "1:10-22:2-14.2:25-75:15-65"; # hsps are sorted by e-value, not position
  #note that coords might be reversed for DNA alignments on reverse strand
  foreach (@hsp)
  {
    if (/(\d+)\:(\d+)\-(\d+)\:(\d+)\-(\d+)/)
    {
      if($5>$4)
      {
        $sub_start{$1}=$4;
        $sub_length{$1}=$5-$4+1;

        if(!$sub_segment_first || $4 < $sub_segment_first){ $sub_segment_first = $4; }
        if(!$sub_segment_last  || $5 > $sub_segment_last){ $sub_segment_last = $5; }
      }
      else # subject in reverse strand
      {
        $sub_start{$1}=$5;
        $sub_length{$1}=$4-$5+1;

        if(!$sub_segment_first || $5 < $sub_segment_first){ $sub_segment_first = $5; }
        if(!$sub_segment_last  || $4 > $sub_segment_last){ $sub_segment_last = $4; }
      }

      $query_start{$1}=$2;
      $query_length{$1}=$3-$2+1;

      if(!$query_segment_first || $2 < $query_segment_first){ $query_segment_first = $2; }
      if(!$query_segment_last  || $3 > $query_segment_last){ $query_segment_last = $3; }
    }
  }

  if($use_short_sequence)
  {
    if($lengthq > $lengths)
    {
      return ( 100*matchlen(\%sub_start,\%sub_length)/$lengths );
    }
    else
    {
      return ( 100*matchlen(\%query_start,\%query_length)/$lengthq );
    }
  }
  else # default get_homologues.pl, calculate cover with respect to longest sequence
  {
    if($lengths > $lengthq)
    {

      #$max_segment = 100* ($sub_segment_last - $sub_segment_first + 1)/$lengths;
      return ( 100*matchlen(\%sub_start,\%sub_length)/$lengths );
    }
    else
    {

       #$max_segment = 100* ($query_segment_last - $query_segment_first + 1)/$lengthq;
      return ( 100*matchlen(\%query_start,\%query_length)/$lengthq );
    }
  }
}

# This subroutine, together with simspan, are used to calculate
# how much of the query sequences match each other.
# Two Arguments:
# 1. Reference Variable: reference to a hash which stores the starting position of each HSP.
# 2. Reference Variable: reference to a hash which stores the length of each HSP.
# Modified by Bruno after INPARANOID: Aug2008,Jan2015
sub matchlen
{
  my ($ref_start,$ref_length) = @_; #10,25 & 13,51

  my @starts = sort{$ref_start->{$a}<=>$ref_start->{$b}} (keys %$ref_start); # (1,2)
  return $ref_length->{$starts[0]} if(scalar(@starts)==1); # 10

  my $i=1;
  my $match_length = $ref_length->{$starts[0]}; # 13
  my $pos = $match_length + $ref_start->{$starts[0]}; # 10

  while($i<scalar(@starts))
  {
    if($ref_length->{$starts[$i]} + $ref_start->{$starts[$i]} <= $pos) # if(51 + 25 <= 10) ????
    {
      $i++;
      next;
    }

    if($ref_start->{$starts[$i]}> $pos){ $match_length += $ref_length->{$starts[$i]}; } # if(25> 10){ match_length = 13 + 51 }
    else{ $match_length += $ref_length->{$starts[$i]} - ($pos - $ref_start->{$starts[$i]}); } # match_length = 13 + 51 - (10 - 25)

    $pos = $ref_start->{$starts[$i]} + $ref_length->{$starts[$i]}; # 25 + 51
    $i++;
  }
  return $match_length;
} ## matchlen

# returns: 1) number of sequences in the array have similar length
# returns: 2) hash with short sequences IDs
# compares length of cluster sequences identified in $ref_array_sequenceIDs with
# respect to first, reference sequence
# uses globals: $MIN_PERCENT_LENGTH_DIFFERENCE if $perc_length not set
# created July 2007, changed July2012
sub same_length
{
  my ($ref_hash_length,$ref_array_sequenceIDs,$perc_length) = @_;

  my ($id,$length,$ref_length,$diff,%short_sequence);
  my $n_of_similar = 1;

  if(!$perc_length){ $perc_length = $MIN_PERCENT_LENGTH_DIFFERENCE }

  $ref_length = $ref_hash_length->{ $ref_array_sequenceIDs->[0] };
  for $id (1 .. $#{$ref_array_sequenceIDs})
  {
    $length = $ref_hash_length->{ $ref_array_sequenceIDs->[$id] };

    if($length > $ref_length){ $diff = $ref_length/$length }
    else{ $diff = $length/$ref_length } #print "#> $id $length $ref_length $diff\n";

    if((100*$diff) >= $perc_length){ $n_of_similar++; }
    else{ $short_sequence{ $ref_array_sequenceIDs->[$id] } = 1 }
  }  #print "$n_of_similar\n";

#	my $bpoline = (split(/;/,$blastquery{$sequenceIDs[0]}))[0];
#	$ref_length = (&getline_from_bpofile($bpoline))[2];
#
#	for($id=1;$id<scalar(@sequenceIDs);$id++)
#	{
#		$bpoline = (split(/;/,$blastquery{$sequenceIDs[$id]}))[0];
#		$length = (&getline_from_bpofile($bpoline))[2];
#
#		if($length > $ref_length){ $diff = $ref_length/$length }
#		else{ $diff = $length/$ref_length } #print "#> $id $length $ref_length $diff\n";
#
#		if((100*$diff) >= $MIN_PERCENT_LENGTH_DIFFERENCE){ $n_of_similar++; }
#		else{ $short_sequence{$sequenceIDs[$id]} = 1 }
#	}

  return ($n_of_similar,\%short_sequence);
} ## same_length

# requires previous call to blast_parse_COG, which will put necessary data in $TMP_DIR
# Updated Aug2016
sub find_COGs
{
  my ($coverage_cutoff,$evalue_cutoff,$multicluster,$saveRAM,$pwd,$force_parsing,@taxa_used) = @_;
  my ($rerun,$clusterid,$id,$genome,$i,$j,$c,$lsematch,$n_of_members,$n_of_clusters) = (0);
  my (%COGs,%orths,%orth_taxa,@LSEclusters,%LSEmembers,%LSEgenome,%singletons,%valid_taxa);

  $coverage_cutoff = sprintf("%1.2f",$coverage_cutoff/100); # re-scale
  if($coverage_cutoff > 0.99)
  { 
    $coverage_cutoff = 0.99; 
    print "\n# find_COGs: coverage limited to -C 99\n\n" ;
  } 

  foreach $genome (@taxa_used){ $valid_taxa{$genome}=1 }

  # check previous runs and save parameters of current run
  if($force_parsing || !-s $parameter_COGS_log || !-e "$TMP_DIR/cog-edges.txt" ||
    check_different_params('COGS',('COVER'=>$coverage_cutoff,
        'EVALUE'=>$evalue_cutoff,'MULTICLUSTER'=>$multicluster))){ $rerun = 1; }

  if($rerun)
  {
    print("# checking lineage-specific expansions\n");
    open (LSEIN,">$lsefilename") || die "# ERROR: find_COGs cannot create $lsefilename, exit\n";
    foreach $i (0 .. $#taxa_used)
    {
      foreach $j (0 .. $#taxa_used)
      {
        next if($i == $j);
        print LSEIN "$taxa_used[$i],$taxa_used[$j]\n";
      }
    }
    close LSEIN;

    # produces _hits, _fhits,_hash.tmp, all.lse.csv
    system("$LSEEXE -t=$TMP_DIR -d=$TMP_DIR -j=$lsefilename -p=$p2ofilename -o=$lseoutfilename > /dev/null");
    if($? != 0)
    {
      die "# ERROR: find_COGs ($LSEEXE) failed to terminate job\n";
      #$LSEEXE -t=$TMP_DIR -d=$TMP_DIR -j=$lsefilename -p=$p2ofilename -o=$lseoutfilename)\n";
    }

    # produces COGtriangles.log,cog-edges.txt,all-edges.txt,all.cog.clusters.log
    print("# making COGs\n");
    
    #print ("$COGTRIANGLESEXE -i=$TMP_DIR -q=$p2ofilename -l=$lseoutfilename -o=$coglogfilename " .
    #        "-n='cluster' -t=$coverage_cutoff -e=$evalue_cutoff\n");
    
    # sleep is to compensate latency/buffering; 
    # NOTE: do not use &> or perl wont wait for system child proc
    system("$COGTRIANGLESEXE -i=$TMP_DIR -q=$p2ofilename -l=$lseoutfilename -o=$coglogfilename " .
        "-n='cluster' -t=$coverage_cutoff -e=$evalue_cutoff 1>/dev/null 2>/dev/null"); sleep(10); 
        
# commented as it takes too much RAM
#if($multicluster){
#	system("$COGTRIANGLESEXE -r -n -i=$TMP_DIR -q=$p2ofilename -l=$lseoutfilename -o=$coglogfilename " .
#       		"-n='cluster' -t=$coverage_cutoff -e=$evalue_cutoff &> /dev/null") }
#else{
#	system("$COGTRIANGLESEXE -r -s -i=$TMP_DIR -q=$p2ofilename -l=$lseoutfilename -o=$coglogfilename " .
#                "-n='cluster' -t=$coverage_cutoff -e=$evalue_cutoff &> /dev/null") }

    # it puts in current directory (pwd) two files (all-edges.txt , cog-edges.txt)
    # Aug2016: confirm that all-edges.txt cannot be empty
    if(-e $pwd.'/cog-edges.txt' && -s $pwd.'/all-edges.txt' && -s $coglogfilename && $? == 0)
    {
      # cog-edges.txt might be actually empty if no COGs are found
      system("mv -f $pwd/all-edges.txt $pwd/cog-edges.txt $TMP_DIR");
    }
    else{ die "# ERROR: find_COGs ($COGTRIANGLESEXE) failed to terminate job\n" }

    save_params('COGS',('COVER'=>$coverage_cutoff,'EVALUE'=>$evalue_cutoff,'MULTICLUSTER'=>$multicluster));
  }

  # parse lineage-specific expansion clusters
  open(LSE,$lseoutfilename) || die "# $0 : cannot read $lseoutfilename\n";
  while(<LSE>)
  {
    #174,175,...
    next if not(/,/);
    push(@LSEclusters,(split)[0]);
  }
  close(LSE);

  # group sequence IDs and count number of taxa in each cluster
  open(COGS,$coglogfilename) || die "# $0 : cannot read $coglogfilename\n";
  while(<COGS>)
  {
    #2,A.faa,2,510,1,510,cluster00001,
    #1,A.faa,1,250,1,250,,
    my @data = split(/,/,$_);
    ($id,$genome,$clusterid) = @data[0,1,6];

    next if(!$valid_taxa{$genome}); # make sure excluded genomes are not considered

    if($clusterid eq '') # singlets in COGtriangles.reformat.pl Jul2012
    {
      $clusterid = ++$n_of_clusters; #print "$clusterid\n";
      $singletons{$id} = $clusterid;

      $lsematch = first { /(?:^|,)$id(?:,|$)/ } @LSEclusters;
      if(defined($lsematch))
      {
        # save observed LSEcluster members
        push(@{$LSEmembers{$lsematch}},$id); #print "$lsematch\n";
        $LSEgenome{$lsematch} = $genome;
      }
    }
    else
    {
      $clusterid = substr($clusterid,7);
      $n_of_clusters = $clusterid;
    }

    #$clusterid =~ /clusterS/ ||  # singletons, deprecated in recent versions of COGtriangles, Jun2012
    #$clusterid =~ /clusterT/); # twogs, frequently contain sequences already in cogs
    push(@{$COGs{$clusterid}{'genomes'}},$genome);
    push(@{$COGs{$clusterid}{'ids'}},$id);
  }
  close(COGS);

  # enforce clusters containing singleton LSEs, otherwise they would be reported as singletons Dic2012
  foreach $c (keys(%LSEmembers))
  {
    $n_of_members = scalar(@{$LSEmembers{$c}});
    next if($n_of_members<2);

    # delete previous singleton clusters
    foreach $id (@{$LSEmembers{$c}})
    {
      delete($COGs{$singletons{$id}});
      #print "deleting singleton $id\n";
    }

    # create new single-genome LSE cluster
    $clusterid = ++$n_of_clusters;
    push(@{$COGs{$clusterid}{'ids'}},@{$LSEmembers{$c}});
    push(@{$COGs{$clusterid}{'genomes'}}, ($LSEgenome{$c}) x $n_of_members);
  }
  undef(@LSEclusters); # free unwanted RAM
  undef(%LSEmembers);
  undef(%LSEgenome);
  undef(%singletons);

  # if required parse clusters and identify sequences included in 1+ clusters
  # just running COGtriangles.reformat.pl wastes too much RAM
  if(!$multicluster) # -s strict mode, sequences assigned to largest clusters
  {
    my (%best_cluster,%best_genome,%best_size,%best_ambiguous);
    my (%multicluster,$genomes,$size);

    # to be used by BerkeleyDB
    my $COG_size_filename  = $TMP_DIR."/COG_size.db";
    my $COG_genome_filename = $TMP_DIR."/COG_genome.db";
    my $COG_cluster_filename  = $TMP_DIR."/COG_cluster.db";
    my $COG_ambiguous_filename = $TMP_DIR."/COG_ambiguous.db";
    my $COG_multicluster_filename = $TMP_DIR."/COG_multicluster.db";

    print("# prunning COGs\n");

    if($saveRAM)
    {
      eval
      {
        import DB_File;
        tie(%best_size,'DB_File',$COG_size_filename,1,0666,$DB_File::DB_HASH)
          || die "# EXIT : cannot create $COG_size_filename: $! (BerkeleyDB::Error)\n";
        tie(%best_genome,'DB_File',$COG_genome_filename,1,0666,$DB_File::DB_HASH)
          || die "# EXIT : cannot create $COG_genome_filename: $! (BerkeleyDB::Error)\n";
        tie(%best_cluster,'DB_File',$COG_cluster_filename,1,0666,$DB_File::DB_HASH)
          || die "# EXIT : cannot create $COG_cluster_filename: $! (BerkeleyDB::Error)\n";
        tie(%best_ambiguous,'DB_File',$COG_ambiguous_filename,1,0666,$DB_File::DB_HASH)
          || die "# EXIT : cannot create $COG_ambiguous_filename: $! (BerkeleyDB::Error)\n";
        tie(%multicluster,'DB_File',$COG_multicluster_filename,1,0666,$DB_File::DB_HASH)
          || die "# EXIT : cannot create $COG_multicluster_filename: $! (BerkeleyDB::Error)\n";
      }
    }

    # find out the largest cluster for each sequence
    foreach $clusterid (keys(%COGs))
    {
      $genomes = scalar(@{$COGs{$clusterid}{'genomes'}});
      $size = scalar(@{$COGs{$clusterid}{'ids'}});
      foreach $id (@{$COGs{$clusterid}{'ids'}})
      {
        $multicluster{$id}++;

        if(!defined($best_cluster{$id}) ||
          $genomes > $best_genome{$id} ||	($genomes == $best_genome{$id} && $size > $best_size{$id}))
        {
          $best_size{$id}      = $size;
          $best_genome{$id}    = $genomes;
          $best_cluster{$id}   = $clusterid;
          $best_ambiguous{$id} = 0;

          #print "> $id $best_size{$id} $best_genome{$id} ".
          #        "$best_cluster{$id} $multicluster{$id} ".
          #	"$best_ambiguous{$id}\n" if($id eq '1001');
        }
        elsif($genomes == $best_genome{$id} && $size == $best_size{$id})
        {
          $best_ambiguous{$id}++;
        }
      }
    }

    # free some hash RAM
    #foreach $id (keys(%best_cluster))
    {
      if($multicluster{$id} < 2)
      {
        delete($best_cluster{$id});
        delete($multicluster{$id});
        delete($best_size{$id});
        delete($best_genome{$id});
        delete($best_ambiguous{$id});
      }

      #else{ print "$id $best_cluster{$id}\n" }
    }

    # make non-redundant clusters
    foreach $clusterid (keys(%COGs))
    {
      my (@nr_ids,@nr_genomes);
      for($i=0;$i<scalar(@{$COGs{$clusterid}{'ids'}});$i++)
      {
        $id = $COGs{$clusterid}{'ids'}[$i];
        if($best_cluster{$id} && ($best_cluster{$id} ne $clusterid || $best_ambiguous{$id}))
        {

          #print "removing sequence $id from cluster $clusterid\n";
          next;
        }
        $genome = $COGs{$clusterid}{'genomes'}[$i];
        push(@nr_ids,$id);
        push(@nr_genomes,$genome);
      } #printf("%s %d %d\n",$clusterid,scalar(@nr_ids),scalar(@nr_genomes));

      if(!@nr_ids)
      {
        delete($COGs{$clusterid});

        #print "deleting cluster $clusterid\n";
      }
      else
      {
        @{$COGs{$clusterid}{'ids'}} = @nr_ids;
        @{$COGs{$clusterid}{'genomes'}} = @nr_genomes;
      }
    }

    if($saveRAM)
    {
      untie(%best_size); untie(%best_genome); untie(%best_cluster);
      untie(%best_ambiguous); untie(%multicluster);
      unlink($COG_size_filename,$COG_genome_filename,$COG_cluster_filename,
        $COG_ambiguous_filename,$COG_multicluster_filename);
    }
    print("# done\n");
  }
  else
  {

    # otherwise do nothing, multicluster is default behaviour
  }

  # finally store all clusters after checking their quality
  CLUSTER: foreach $clusterid (sort {$a<=>$b} (keys(%COGs)))
  {
    my @cluster = sort {$a<=>$b} (@{$COGs{$clusterid}{'ids'}});
    ($id,$i) = (shift(@cluster),1);
    while($orths{$id})
    {
      push(@cluster,$id);
      if($i == scalar(@cluster))
      {
        print "# find_COGs : cannot find unique id for cluster $clusterid, skip it\n";
        next CLUSTER;
      }
      $id = shift(@cluster);
      $i++;
    }
    push(@{$orths{$id}},@cluster);

    foreach $genome (@{$COGs{$clusterid}{'genomes'}})
    {
      $orth_taxa{$id}{$genome}++;
    }
  }

  return(\%orths,\%orth_taxa);
}

# uses globals: %gindex , Oct2015
sub flag_small_clusters
{
  my ($ref_hash_orths,$ref_hash_orth_taxa,$ref_used_taxa,$min_cluster_size) = @_;

  my ($id,$first,$last,%small,%large,%stats);

  # store sequence ids of clusters >= than required
  foreach $id (keys(%$ref_hash_orths))
  {
    next if(scalar(keys(%{$ref_hash_orth_taxa->{$id}})) < $min_cluster_size);
    $large{$id}=1;
    foreach my $m (@{$ref_hash_orths->{$id}})
    {
      $large{$m}=1;
    }
  }

  # loop across sequence ids in each taxon and flag those in small clusters
  foreach my $taxon (@$ref_used_taxa)
  {
    $first = $gindex{$taxon}[0];
    $last  = $gindex{$taxon}[1];
    foreach $id ($first .. $last)
    {
      next if($large{$id});
      $small{$id} = 1;
      $stats{$taxon}++;
    }
  }

  print "\n# flag_small_clusters:\n";
  foreach my $taxon (@$ref_used_taxa)
  {
    print "# $taxon : $stats{$taxon}\n";
  }
  printf("# total : %d sequences\n\n",scalar(keys(%small)));

  return \%small;
}


# adds to orthologous clusters single sequences with no valid BLAST hits, so that
# pangenomes are complete Zhang,Yunzeng
# created by Bruno mar2014, fixed nov2014, updated Oct2015
# uses globals: %gindex,@gindex2
sub add_unmatched_singletons
{
  my ($ref_hash_orths,$ref_hash_orth_taxa,$ref_used_taxa,$ref_sequences2skip) = @_;

  my ($id,$prev,$first,$last,$taxon);
  my (@unmatched,%added_taxa);
  my $first_id = 9999999999999999999;
  my $last_id = 0;

  # find out range of valid sequence ids
  foreach $taxon (@$ref_used_taxa)
  {
    $first = $gindex{$taxon}[0];
    $last  = $gindex{$taxon}[1];
    if($first < $first_id){ $first_id = $first }
    if($last > $last_id ){ $last_id = $last }
  }

  # sort sequence ids included in clusters and find gaps (unmatched ids)
  my ($tmpfh,$tmpfilename) = tempfile(UNLINK=>1);
  foreach $id (keys(%$ref_hash_orths))
  {
    print $tmpfh "$id\n";
    foreach my $m (@{$ref_hash_orths->{$id}})
    {
      print $tmpfh "$m\n";
    }
  }
  close($tmpfh);

  $prev = $first_id - 1; # fixed Jan2015: -1
  $ENV{'LC_ALL'} = 'POSIX';
  open(SORT,"$SORTBIN -k 1,1n $tmpfilename |")
    || die "# add_unmatched_singletons : cannot sort\n";
  while(<SORT>)
  {
    $id = (split)[0];
    while(($id-$prev)>1){ $prev++; push(@unmatched,$prev); }

    $prev = $id;
  }
  close(SORT);

  while($id < $last_id)
  {
    $id++;
    push(@unmatched,$id)
  }

  # add unmatched ids to input data structures
  my $actually_added = 0;
  foreach $id (@unmatched)
  {
    $taxon = $gindex2[$id];
    next if(!grep(/^$taxon$/,@$ref_used_taxa));
    next if($ref_sequences2skip->{$id});

    $added_taxa{$taxon}++; #print "$id $taxon\n";

    push(@{$ref_hash_orths->{$id}},());
    $ref_hash_orth_taxa->{$id}{$taxon}++;
    $actually_added++;
  }

  print "\n# add_unmatched_singletons : $actually_added sequences, ".
    scalar(keys(%added_taxa))." taxa\n\n"; #print "$tmpfilename\n";
}

# Runs PARANOID algorithm to cluster orthologs and their inparalogs
# adapted by Bruno Aug2008
sub find_PARANOID_clusters
{
  my ($saveRAM,$bitscore_cutoff,$pmatch_cutoff,$psegment_cutoff,$path,$force_parsing,@taxa) = @_;

  my ($i,$j,@inparanoids,$inparaoutfile);
  my (%raw_clusters,%orthologues,%orth_taxa,$cluster,$orth,$ref,$taxon);
  my $rerun = 0;
  my $multiparanoid_command = "$MULTIPARANOID ";

  # 0) check previous runs and save parameters of current run
  if($force_parsing || !-s $parameter_PRND_log ||
    check_different_params('PRND',('B'=>$bitscore_cutoff,'COVER'=>$pmatch_cutoff,'SEG'=>$psegment_cutoff)))
  {
    $rerun = 1;
  }

  save_params('PRND',('B'=>$bitscore_cutoff,'COVER'=>$pmatch_cutoff,'SEG'=>$psegment_cutoff));

  # 1) run inparanoid with all possible pairs of taxa
  for($i=0;$i<scalar(@taxa)-1;$i++)
  {
    for($j=$i+1;$j<scalar(@taxa);$j++)
    {
      print("\n# Identifying INPARANOID orthologous pairs between $i:$taxa[$i] and $j:$taxa[$j] ...\n");
      $inparaoutfile =
        execute_INPARANOID($saveRAM,$taxa[$i],$taxa[$j],$i,$j,$bitscore_cutoff,$pmatch_cutoff,$psegment_cutoff,$path,$rerun);
      push(@inparanoids,$inparaoutfile)
    }
  }

  # 2) run multiparanoid to merge clusters
  if($rerun || !-s $paranoid_file)
  {
    my $outmultiparanoidfile = '';

    for($i=0;$i<scalar(@taxa);$i++){ $multiparanoid_command .= "$i+"; }
    chop($multiparanoid_command);

    open(MULTI,"$multiparanoid_command $path 2>&1 |"); #print ">> $multiparanoid_command $path\n";
    while(<MULTI>)
    {
      if(/The result is in file (\S+)/){ $outmultiparanoidfile = $1 } #print;
    }
    close(MULTI);

    if(!$outmultiparanoidfile)
    {
      die "# find_PARANOID_clusters : could not find parse multiparanoid outfile! ($multiparanoid_command $path)\n";
    }
    else
    {
      system("cp $outmultiparanoidfile $paranoid_file");
      print("\n# multiparanoid: $multiparanoid_command $path\n\n");
    }
  }
  else{ print("\n# re-using previous multiparanoid results...\n\n"); }

  # 3) parse multiparanoid outfile
  open(MULTIOUT,$paranoid_file) ||
    die "# find_PARANOID_clusters : cannot find $paranoid_file, probably a multiparanoid error! ($multiparanoid_command $path)\n";
  while(<MULTIOUT>)
  {
    #clusterID	species	gene	is_seed_ortholog	confidence_score	species_in_cluster	tree_conflict
    #1	Ath.faa	6	1	1.000	Ath.faa-Hsa.faa	No
    if(/(\d+)\t(\S+)\t(\S+)/)
    {
      push(@{$raw_clusters{$1}},$3); #print "#$1 $2 $3\n";
    }
  }
  close(MULTIOUT);

  # 4) pack clusters and taxa
  foreach $cluster (keys(%raw_clusters))
  {
    $ref = $raw_clusters{$cluster}[0];

    foreach $orth (@{$raw_clusters{$cluster}})
    {
      $taxon = $gindex2[$orth]; #print "|$ref $orth $taxon|\n";
      $orth_taxa{$ref}{$taxon}++;
    }

    shift(@{$raw_clusters{$cluster}});
    foreach $orth (@{$raw_clusters{$cluster}}) { push(@{$orthologues{$ref}},$orth); }
  }

  return (\%orthologues,\%orth_taxa);

} ## find_PARANOID_clusters

# runs INPARANOID algorithm with a pair A,B of taxa
# the code is an adaptation of inparanoid 3.0
# returns: sql table for multiparanoid use
sub execute_INPARANOID
{

  # INPARANOID adapted code removed, authors do not license redistribution
}

# Split clusters of orths generated by OMCL, COG or PRND according to Pfam domain composition
# Updates $ref_orth_taxa and returns \%orthologues
# Uses globals: %pfam_hash, @gindex2
sub split_Pfam_clusters
{
  my ($ref_hash_orths, $ref_orth_taxa, $verbose) = @_;

  my ($pfam,$orth,$subcluster,$cluster,$taxon,%orthologues);
  my ($n_of_original_clusters,$n_of_new_clusters) = (0,0);

  print "# splitting clusters by Pfam domain composition\n";

  foreach $cluster (keys(%$ref_hash_orths))
  {

    # check number of Pfam domain strings in this cluster
    my (%Pfam_strings);

    if(defined($pfam_hash{$cluster})){ $pfam = "$pfam_hash{$cluster}" }
    else{ $pfam = '' }
    push(@{$Pfam_strings{$pfam}},$cluster);
    foreach $orth (@{$ref_hash_orths->{$cluster}})
    {
      if(defined($pfam_hash{$orth})){ $pfam = "$pfam_hash{$orth}" }
      else{ $pfam = '' }
      push(@{$Pfam_strings{$pfam}},$orth);
    } #print "<".$#{$ref_hash_orths->{$cluster}}." ".scalar(@{$ref_hash_orths->{$cluster}})."\n";

    # create as many subclusters as diff Pfam strings found
	if(scalar(keys(%Pfam_strings)) > 1)
    {
      $n_of_original_clusters++;

      print "# old cluster $cluster contains these Pfam strings: ".
        join(' ; ',keys(%Pfam_strings))."\n" if($verbose);

      foreach $pfam (keys(%Pfam_strings))
      {
        $n_of_new_clusters++;

        $subcluster = shift(@{$Pfam_strings{$pfam}});
        $taxon = $gindex2[$subcluster];

        @{$orthologues{$subcluster}} = ();
        undef(%{$ref_orth_taxa->{$subcluster}});
        $ref_orth_taxa->{$subcluster}{$taxon}++;

        foreach $orth (@{$Pfam_strings{$pfam}})
        {
          $taxon = $gindex2[$orth];
          push(@{$orthologues{$subcluster}},$orth);
          $ref_orth_taxa->{$subcluster}{$taxon}++;
        } #print ">> $subcluster ".scalar(@{$orthologues{$subcluster}})."\n";
      }
    }
    else # add cluster as is if only one domain combination
    {
      $orthologues{$cluster} = $ref_hash_orths->{$cluster};
    }
  }

  printf("# split_Pfam_clusters : created %d new clusters\n",$n_of_new_clusters-$n_of_original_clusters);

  return \%orthologues;
}

# runs all step required for OMCL clustering trying to re-use as much as possible previous results
# created by Bruno oct2011, updated Jan2015
sub find_OMCL_clusters
{
  my ( $saveRAM,$evalue_cutoff,$pi_cutoff,$pmatch_cutoff,$neighbor_corr_cutoff,
    $MCLinflation,$sorted_inparalogues,$redo_inp,$redo_orth,$redo_mcl,
    $ref_taxa,$ref_full_sequence_taxa ) = @_;

  # to be used by BerkeleyDB
  my $MCL_graph_filename  = $TMP_DIR."/MCL_graph.db";
  my $MCL_weight_filename = $TMP_DIR."/MCL_weight.db";

  if(!$redo_inp && !$redo_orth && !$redo_mcl && -s $mcl_file)
  {
    print "# find_OMCL_clusters: re-using previous results (same parameters)\n";
  }
  else
  {

    # rewrite matrx file if any of these is true
    if($redo_inp || $redo_orth || $redo_mcl || !-s $matrix_file)
    {
      if($saveRAM)
      {
        eval
        {
          import DB_File;
          tie(@graph,'DB_File',$MCL_graph_filename,1,0666,$DB_File::DB_RECNO)
            || die "# EXIT : cannot create $MCL_graph_filename: $! (BerkeleyDB::Error)\n";
          tie(%weight,'DB_File',$MCL_weight_filename,1,0666,$DB_File::DB_HASH)
            || die "# EXIT : cannot create $MCL_weight_filename: $! (BerkeleyDB::Error)\n";

          }
      }

      if($redo_inp == -1){ $redo_inp = 0 }
      if($redo_orth == -1){ $redo_orth = 0 }

      if($ref_full_sequence_taxa)
      {
        findAllOrthologiesORTHMCL($saveRAM,$evalue_cutoff,$pi_cutoff,$pmatch_cutoff,$neighbor_corr_cutoff,
          $sorted_inparalogues,$redo_inp,$redo_orth,$ref_taxa,$ref_full_sequence_taxa);
      }
      else
      {
        findAllOrthologiesORTHMCL($saveRAM,$evalue_cutoff,$pi_cutoff,$pmatch_cutoff,$neighbor_corr_cutoff,
          $sorted_inparalogues,$redo_inp,$redo_orth,$ref_taxa);
      }

      write_matrix_index();

      if($saveRAM)
      {
        untie (@graph); untie(%weight);
        unlink ($MCL_graph_filename,$MCL_weight_filename);
      }
      else{ undef(@graph); undef(%weight); }

      executeMCL($matrix_file,$mcl_file,$MCLinflation);
    }
    elsif(!-s $mcl_file)
    {
      executeMCL($matrix_file,$mcl_file,$MCLinflation);
    }
  }

  print "\n# find_OMCL_clusters: parsing clusters ($mcl_file)\n";
  my ($ref_hash_orths,$ref_hash_orth_taxa) = parse_MCL_clusters( $mcl_file );
  print "\n";

  return ($ref_hash_orths,$ref_hash_orth_taxa);
}

# returns: adds all orthologies found by orth all against all comparisons
# uses globals: %graph,%weight,$bpofile,$last_graph_item
# modified from original orthomcl.pl, last change in May2012
# jan2015 added $ref_full_sequence_taxa
sub findAllOrthologiesORTHMCL
{
  my ($saveRAM, $pv_cutoff, $pi_cutoff, $pmatch_cutoff, $neighbor_corr_cutoff,
    $inparalogues_are_sorted, $redo_inparal, $redo_orth, $ref_taxa, $ref_full_sequence_taxa) = @_;

  my ($p,$i,$j,$n,$k,$partial_sequences,@ortho,@taxa);
  my ($avgw,$ref,$refi);

  @taxa = @$ref_taxa;

  for($i=0;$i<scalar(@taxa);$i++)
  {
    my (@MCL_nodes,@MCL_innodes,%MCL_scores,%scores);
    for($j=$i+1;$j<scalar(@taxa);$j++)
    {
      $partial_sequences = 0;

      if($ref_full_sequence_taxa &&
        (!$ref_full_sequence_taxa->{$taxa[$i]} || !$ref_full_sequence_taxa->{$taxa[$j]})){ $partial_sequences = 1 }

      print("\n# identifying orthologs between $taxa[$i] and $taxa[$j] ($partial_sequences)\n");

      @{$ref} = makeOrtholog($saveRAM,$taxa[$i],$taxa[$j],0,$pv_cutoff,$pi_cutoff,
        $pmatch_cutoff,$neighbor_corr_cutoff,1,0,$redo_inparal || $redo_orth, {}, {}, $partial_sequences);

      my $edge_ref   = $ref->[0];
      my %w          = %{$ref->[1]};
      my $sumw       = $ref->[3];
      my $c_ortholog = $ref->[4];
      my @nodes = sort {$a<=>$b} keys %{$edge_ref};
      foreach $n (@nodes)
      {
        $ortho[$n] = 1;
        $MCL_nodes[$n] = join(',',@{$edge_ref->{$n}});
      }
      if($nodes[$#nodes] > $last_graph_item){ $last_graph_item = $nodes[$#nodes] }

      # save normalized weights
      if($c_ortholog>0)
      {
        $avgw = $sumw/$c_ortholog;
        foreach $p (keys(%w)){ $MCL_scores{$p} = sprintf("%.3f", $w{$p}/$avgw); }
      }
      else{ foreach $p (keys(%w)){ $MCL_scores{$p} = $w{$p} } }
    }

    # normalize inparalogue's weights
    my ($count,$sum,$count_all,$sum_all) = (0,0,0);
    print("\n# identifying inparalogs in $taxa[$i]\n");
    @{$refi} = makeInparalog($saveRAM,$taxa[$i],$pv_cutoff,$pi_cutoff,$pmatch_cutoff,
      $neighbor_corr_cutoff,$inparalogues_are_sorted,$redo_inparal,$partial_sequences);
    my @innodes = sort {$a<=>$b} keys %{$refi->[0]};
    foreach $p (@innodes) { $MCL_innodes[$p] = join(',',@{$refi->[0]->{$p}}); }
    if($innodes[$#innodes] > $last_graph_item){ $last_graph_item = $innodes[$#innodes] }
    foreach $p (keys %{$refi->[1]}) { $scores{$p} = $refi->[1]->{$p}; } # 2 concatenated keys
    my @pairs = keys(%scores);
    foreach $p (@pairs)
    {
      ($n,$k) = split(' ',$p);
      $count_all++;
      $sum_all += $scores{$p};
      if ($ortho[$n] || $ortho[$k])
      {
        $count++;
        $sum += $scores{$p};
      }
    }

    # normalize inparalog weights by the average weight of inparalogs which have orthologs in other species
    # common case, for eukaryotes and most prokaryotes
    if($count)
    {
      $avgw = $sum/$count;
      foreach $p (@pairs){ $MCL_scores{$p} = sprintf("%.3f", $scores{$p}/$avgw) }
    }

    # OR normalize the in-paralog weights by the average weight of all inparalogs
    # not common, useful for prokaryotes or pathogens
    elsif($count_all)
    {
      $avgw = $sum_all/$count_all;
      foreach $p (@pairs){ $MCL_scores{$p} = sprintf("%.3f", $scores{$p}/$avgw) }
    }

    # OR no normalization since $count_all=0 and there is nothing stored in %weight
    # not common, useful for prokaryotes or pathogens
    else { foreach $p (@pairs){ $MCL_scores{$p} = sprintf("%.3f", $scores{$p})  }}

    # put together calculated nodes and scores in globals
    foreach $p (keys(%MCL_scores))
    {
      $weight{$p} = $MCL_scores{$p}; #printf("%s %f\n",$p,$MCL_scores{$p});
    }
    foreach $p (0 .. $#MCL_nodes)
    {
      next if(!$MCL_nodes[$p]); #printf(":%d:  %s,\n",$p,$MCL_nodes[$p]);
      $graph[$p] .= "$MCL_nodes[$p],";
    }
    foreach $p (0 .. $#MCL_innodes)
    {
      next if(!$MCL_innodes[$p]); #printf(":%d:  %s,\n",$p,$MCL_innodes[$p]);
      $graph[$p] .= "$MCL_innodes[$p],";
    }
  }
} ## findAllOrthologiesORTHMCL

# returns: [-1,1] correlation, inspired by Song,Joseph,Davis,Durand(2008)PlOS PMID:18475320
# uses globals: %blastquery (caution if $saveRAM is used!)
# takes log10(bit-scores) for calculations
# decimal differences have large impact on correlation calculations with small samples, check printf("%1.0f",...)
# Created by Bruno Aug2008
sub neighbor_correlation
{
  my ($seqid1,$seqid2) = @_; #print "|$seqid1|$seqid2|\n";

  my ($NC,$num,$denom,$var1,$var2,$mean1,$mean2,$total) = (0,0,0,0,0,0,0,0);
  my ($bits,$line,$sid,%bitsim) = (-1,-1);
  my $unrelated_bitsim =  sprintf("%1.0f", log($MIN_BITSCORE_SIM_DEFAULT)/log(10) );

# 1) recorrer lista de hits de 1, guardar vector de similitudes y calcular similitud media
  if($blastquery{$seqid1})
  {
    my $ref_lines = getblock_from_bpofile($blastquery{$seqid1}[0],$blastquery{$seqid1}[1]);
    foreach $line (@{$ref_lines})
    {
      #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      ($sid,$bits) = @$line[1,9];
      $bitsim{$sid}{$seqid1} =  sprintf("%1.0f", log($bits)/log(10) );
    }
  }

  # 2) recorrer lista de hits de 2, guardar vector de similitudes y calcular similitud media
  if($blastquery{$seqid2})
  {
    my $ref_lines = getblock_from_bpofile($blastquery{$seqid2}[0],$blastquery{$seqid2}[1]);
    foreach $line (@{$ref_lines})
    {
      ($sid,$bits) = @$line[1,9];
      $bitsim{$sid}{$seqid2} =  sprintf("%1.0f", log($bits)/log(10) ); #print "$seqid2|$sid|$bits|$bitsim{$sid}{$seqid2}|\n";
    }
  }

# 3) calcular las similitudes medias de ambas secuencias, incluyendo las secuencias no relacionadas
  foreach $sid (sort {$a<=>$b} (keys(%bitsim)))
  {
    if(!$bitsim{$sid}{$seqid1}){ $bitsim{$sid}{$seqid1} = $unrelated_bitsim; }
    elsif(!$bitsim{$sid}{$seqid2}){ $bitsim{$sid}{$seqid2} = $unrelated_bitsim; } #print "$sid\t$bitsim{$sid}{$seqid1}\t$bitsim{$sid}{$seqid2}\n";
    $mean1 += $bitsim{$sid}{$seqid1};
    $mean2 += $bitsim{$sid}{$seqid2};
    $total++;
  }
  $mean1 /= $total;
  $mean2 /= $total; #print "#$total#$mean1#$mean2#\n";

  # 4) calcular numerador de ecuacion1 del paper, teniendo en cuenta N y similitud por defecto
  foreach $sid (keys(%bitsim)){ $num += (($bitsim{$sid}{$seqid1} - $mean1) * ($bitsim{$sid}{$seqid2} - $mean2)); }

  # 5) calcular varianzas de similitud de 1 y 2
  foreach $sid (keys(%bitsim))
  {
    $var1 += ($bitsim{$sid}{$seqid1} - $mean1) ** 2;
    $var2 += ($bitsim{$sid}{$seqid2} - $mean2) ** 2;
  }

  # 7) calcular denominador de ecuacion1
  $denom = sqrt( $var1 * $var2 );

  # 8) calcular neighbor correlation (NC)
  if($denom){ $NC = sprintf("%1.2f",$num/$denom); }
  elsif($num == 0 && $denom == 0){ $NC = 1.00 }        #print "#$num#$denom#$NC#\n";

  return $NC;
} ## neighbor_correlation

# returns a hash with first and last residues aligned in observed local alignments
# Bruno Oct2015
sub find_local_alignment_coords
{
  my (@ids) = @_;

  my ($id,$sid,$line,$simspan,$length,$first,$last,%valid_subjects,%coords);

  foreach $id (@ids){ $valid_subjects{$id} = 1 }

  foreach $id (@ids)
  {
    my $ref_lines = getblock_from_bpofile($blastquery{$id}[0],$blastquery{$id}[1]);
    foreach $line (@{$ref_lines})
    {
      #$pQid\t$pSid\t$pEvalue\t$ppercID\t$Qcov\t$Scov\t$pQlength\t$pSlength\t$simspan\t$pbits
      ($sid,$length,$simspan) = @$line[1,6,8];
      next if($sid == $id || !$valid_subjects{$sid});

      ($first,$last) = alignment_coords( $simspan );

      if(!$coords{$id}{'first'} || $first < $coords{$id}{'first'}){ $coords{$id}{'first'} = $first }
      if(!$coords{$id}{'last'} || $last > $coords{$id}{'last'}){ $coords{$id}{'last'} = $last }
    }

    $coords{$id}{'length'} = $length;
  }

  return %coords;
} ## find_local_alignment_coords

# check genome files used to generate a given all.bpo file ($selected_genomes_file)
sub get_string_with_previous_genomes
{
  my ($selected_genomes_file) = @_;

  my ($previous_genomes,@genomes) = ('');

  open(SEL,$selected_genomes_file) || die "# get_string_with_previous_genomes : cannot read $selected_genomes_file\n";
  while(<SEL>)
  {
    chomp $_;
    $_ =~ s/\s+//g;
    push(@genomes,$_);
  }
  close(SEL);

  $previous_genomes = join('',sort(@genomes));

  return $previous_genomes;
}

1;
