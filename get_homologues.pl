#!/usr/bin/env perl 

# 2017 Bruno Contreras-Moreira (1) and Pablo Vinuesa (2):
# 1: http://www.eead.csic.es/compbio (Estacion Experimental Aula Dei/CSIC/Fundacion ARAID, Spain)
# 2: http://www.ccg.unam.mx/~vinuesa (Center for Genomic Sciences, UNAM, Mexico)

# This program uses BLAST (and DIAMOND/HMMER/Pfam) to define clusters of 'orthologous' ORF/intergenic
# sequences and pan/core-genome gene sets. Several algorithms are available and search parameters
# are customizable. It is designed to process (in a SGE computer cluster) files contained in a
# directory (-d), so that new .faa/.gbk files can be added while conserving previous BLAST results.
# In general the program tries to re-use previous results when run with the same input directory.

$|=1;

use strict;
use Getopt::Std;
use Fcntl qw(:flock);
use File::Temp qw(tempfile);
use File::Basename;
use Benchmark;
use Cwd;
use FindBin '$Bin';
use lib "$Bin/lib";
use lib "$Bin/lib/bioperl-1.5.2_102/";
use phyTools; #also imports constants SEQ,NAME used as array subindices
use marfil_homology; # includes $FASTAEXTENSION and other global variables set there such as @taxa

my $VERSION = '3.x'; 

## sun grid engine (computer cluster) variables, might require edition to fit your system (ignored with -m local)
my $SGEPATH = "";
if($ENV{'SGE_ROOT'} && $ENV{'SGE_ARCH'}){ $SGEPATH = $ENV{'SGE_ROOT'}.'/bin/'.$ENV{'SGE_ARCH'}; };
my $SGEERRORSTATUS = 'Eqw';       # error status of failed jobs in SGE
my $QUEUESETTINGS = '';
my $QUEUEWAIT = 1;                # interval in seconds between qsub submissions
my $WAITTIME  = 30;               # interval in seconds btween qstat commands

## settings for local batch blast,hmmscan jobs
my $BATCHSIZE = 100;

## global variables that control some algorithmic choices
my $NOFSAMPLESREPORT     = 10;    # number of genome samples used for the generation of pan/core genomes
my $MAXEVALUEBLASTSEARCH = 0.01;  # required to reduce size of blast output files
my $MAXPFAMSEQS          = 250;   # default for -m cluster jobs, it is updated to 1000 with -m local
my $MININTERGENESIZE     = 200;   # minimum length (nts) required for intergenic segments to be considered
my $MAXINTERGENESIZE     = 700;
my $INTERGENEFLANKORF    = 180;   # length in nts of intergene flanks borrowed from neighbor ORFs (/3 for amino acids)
my $PRINTCLUSTERSSCREEN  = 0;     # whether cluster names should be printed to screen
my $KEEPSCNDHSPS         = 1;     # by default only the first BLAST hsp is kept for further calculations
my $DISSABLEPRND         = 1;     # dissable paranoid (PRND, -P) options due to no redistribution license;


## list of features/binaries required by this program (do not edit)
my @FEATURES2CHECK = ('EXE_BLASTP','EXE_BLASTN','EXE_FORMATDB','EXE_MCL','EXE_MAKEHASH',
  'EXE_READBLAST','EXE_COGLSE','EXE_COGTRI','EXE_HMMPFAM','EXE_DMNDP','EXE_DMNFT',
  'EXE_INPARA','EXE_ORTHO','EXE_HOMOL','EXE_SPLITBLAST','EXE_SPLITHMMPFAM');
if(!$DISSABLEPRND){ push(@FEATURES2CHECK,'EXE_MULTIPARANOID'); }

################################################################

my ($newDIR,$output_mask,$pancore_mask,$include_file,$onlyblast,%included_input_files,%opts) = ('','','',0,0);
my ($exclude_inparalogues,$doMCL,$doPARANOID,$doCOG,$do_PFAM,$reference_proteome_string) = (0,0,0,0,0,0);
my ($inputDIR,$input_FASTA_file,$filter_by_length,$cluster_list_file,$do_intergenic,$do_features) = (0,0);
my ($n_of_cpus,$COGmulticluster,$do_minimal_BDBHs,$do_ANIb_matrix,$do_soft) = ($BLAST_NOCPU,0,0,0,0);
my ($do_diamond) = (0);
my ($min_cluster_size,$runmode,$do_genome_composition,$saveRAM,$ANIb_matrix_file);
my ($evalue_cutoff,$pi_cutoff,$pmatch_cutoff) = ($BLAST_PVALUE_CUTOFF_DEFAULT, # see marfil_homology.pm
  $PERCENT_IDENTITY_CUTOFF_DEFAULT,$PERCENT_MATCH_CUTOFF_DEFAULT);
my ($MCLinflation,$neighbor_corr_cutoff) = ($MCL_INFLATION_DEFAULT,$MIN_NEIGHBOR_CORR_DEFAULT);
my ($bitscore_cutoff,$segment_cover_cutoff) = ($MIN_BITSCORE_DEFAULT,$MIN_SEGMENT_COVER_DEFAULT);
my $random_number_generator_seed = 0;
my $pwd = getcwd(); $pwd .= '/';
my $start_time = new Benchmark();

getopts('hvbcesgADGPMXoxza:f:n:m:d:r:t:i:I:E:S:C:F:N:B:O:R:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-v print version, credits and checks installation\n";
  print   "-d directory with input FASTA files ( .faa / .fna ),           (overrides -i,\n";
  print   "   GenBank files ( .gbk ), 1 per genome, or a subdirectory      use of pre-clustered sequences\n".  
          "   ( subdir.clusters / subdir_ ) with pre-clustered sequences   ignores -c, -g)\n".
          "   ( .faa / .fna ); allows for new files to be added later;    \n".
          "   creates output folder named 'directory_homologues'\n";
  print   "-i input amino acid FASTA file with [taxon names] in headers,  ".
    "(required unless -d is set)\n".
    "   creates output folder named 'file_homologues'\n";
  print   "\nOptional parameters:\n";
  print   "-o only run BLAST/Pfam searches and exit                       ".
    "(useful to pre-compute searches)\n";
  print   "-c report genome composition analysis                          ".
    "(follows order in -I file if enforced,\n".
    "                                                               ".
    " ignores -r,-t,-e)\n";
  print   "-R set random seed for genome composition analysis             ".
    "(optional, requires -c, example -R 1234,\n".
    "                                                               ".
    " required for mixing -c with -c -a runs)\n";
  if(eval{ require DB_File } ) # show only if DB_File is available (it should as it is core)
  {
    print   "-s save memory by using BerkeleyDB; default parsing stores\n".
      "   sequence hits in RAM\n";
  }
  print   "-m runmode [local|cluster]                                     ".
    "(default local)\n";
  print   "-n nb of threads for BLAST/HMMER/MCL in 'local' runmode        (default=$n_of_cpus)\n";
  print   "-I file with .faa/.gbk files in -d to be included              ".
    "(takes all by default, requires -d)\n";
  print   "\nAlgorithms instead of default bidirectional best-hits (BDBH):\n";
  print   "-G use COGtriangle algorithm (COGS, PubMed=20439257)           ".
    "(requires 3+ genomes|taxa)\n";
  print   "-M use orthoMCL algorithm (OMCL, PubMed=12952885)\n";
  print   "-P use PARANOID algorithm (PRND, PubMed=16873526)\n" if(!$DISSABLEPRND);
  print   "\nOptions that control sequence similarity searches:\n";
  print   "-X use diamond instead of blastp                               (optional, set threads with -n)\n";
  print   "-C min \%coverage in BLAST pairwise alignments                  ";
  print "(range [1-100],default=$PERCENT_MATCH_CUTOFF_DEFAULT)\n";
  print   "-E max E-value                                                 ".
    "(default=$BLAST_PVALUE_CUTOFF_DEFAULT,max=$MAXEVALUEBLASTSEARCH)\n";
  print   "-D require equal Pfam domain composition                       ".
    "(best with -m cluster or -n threads)\n";
  print   "   when defining similarity-based orthology\n";
  print   "-S min \%sequence identity in BLAST query/subj pairs            ".
    "(range [1-100],default=$PERCENT_IDENTITY_CUTOFF_DEFAULT [BDBH|OMCL])\n";
  print   "-N min BLAST neighborhood correlation PubMed=18475320          ".
    "(range [0,1],default=$MIN_NEIGHBOR_CORR_DEFAULT [BDBH|OMCL])\n";
  print   "-b compile core-genome with minimum BLAST searches             ".
    "(ignores -c [BDBH])\n";
  print   "\nOptions that control clustering:\n";
  print   "-t report sequence clusters including at least t taxa          ".
    "(default t=numberOfTaxa,\n".
    "                                                               ".
    " t=0 reports all clusters [OMCL|COGS])\n";
  print   "-a report clusters of sequence features in GenBank files       ".
    "(requires -d and .gbk files,\n";
  print   "   instead of default 'CDS' GenBank features                   ".
    " example -a 'tRNA,rRNA',\n".
    "                                                               ".
    " NOTE: uses blastn instead of blastp,\n".
    "                                                               ".
    " ignores -g,-D)\n";

  print   "-g report clusters of intergenic sequences flanked by ORFs     ".
    "(requires -d and .gbk files)\n";
  print   "   in addition to default 'CDS' clusters\n";
  print   "-f filter by \%length difference within clusters                ".
    "(range [1-100], by default sequence\n".
    "                                                               ".
    " length is not checked)\n";
  print   "-r reference proteome .faa/.gbk file                           ".
    "(by default takes file with\n".
    "                                                               ".
    " least sequences; with BDBH sets\n".
    "                                                               ".
    " first taxa to start adding genes)\n";
  print   "-e exclude clusters with inparalogues                          ".
    "(by default inparalogues are\n".
    "                                                               ".
    " included)\n";
  print   "-x allow sequences in multiple COG clusters                    ".
    "(by default sequences are allocated\n".
    "                                                               ".
    " to single clusters [COGS])\n";
  print   "-F orthoMCL inflation value                                    ".
    "(range [1-5], default=$MCL_INFLATION_DEFAULT [OMCL])\n";
  print   "-B min bit-score (in bits)                                     ".
    "(default=$MIN_BITSCORE_DEFAULT [PRND])\n" if(!$DISSABLEPRND);
  print   "-A calculate average identity of clustered sequences,          ".
    "(optional, creates tab-separated matrix,\n";
  print   " by default uses blastp results but can use blastn with -a     ".
    " recommended with -t 0 [OMCL|COGS])\n";
  print   "-z add soft-core to genome composition analysis                ".
    "(optional, requires -c [OMCL|COGS])\n";
    
  print "\n".
    " This program uses BLAST (and optionally HMMER/Pfam) to define clusters of 'orthologous'\n".
    " genomic sequences and pan/core-genome gene sets. Several algorithms are available\n".
    " and search parameters are customizable. It is designed to process (in a SGE computer\n".
    " cluster) files contained in a directory (-d), so that new .faa/.gbk files can be added\n".
    " while conserving previous BLAST results. In general the program will try to re-use\n".
    " previous results when run with the same input directory.\n\n";
  exit;
}

# read version bumber from CHANGES.txt
open(CHANGES,"$Bin/CHANGES.txt");
while(<CHANGES>)
{
  if(eof && /^(\d+):/){ $VERSION = $1 } 
}
close(CHANGES);

if(defined($opts{'v'}))
{
  print "\n$0 version $VERSION\n";
  print "\nProgram written by Bruno Contreras-Moreira (1) and Pablo Vinuesa (2).\n";
  print "\n 1: http://www.eead.csic.es/compbio (Estacion Experimental Aula Dei/CSIC/Fundacion ARAID, Spain)\n";
  print " 2: http://www.ccg.unam.mx/~vinuesa (Center for Genomic Sciences, UNAM, Mexico)\n";
  print "\nPrimary citations:\n\n";
  print " Contreras-Moreira B, Vinuesa P (2013) GET_HOMOLOGUES, a versatile software package for scalable and\n".
        " robust microbial pangenome analysis. Appl Environ Microbiol 79(24):7696-701. (PubMed:24096415)\n\n";
  print " Vinuesa P and Contreras-Moreira B (2015) Robust Identification of Orthologues and Paralogues for \n".
        " Microbial Pan-Genomics Using GET_HOMOLOGUES: A Case Study of pIncA/C Plasmids. In Bacterial Pangenomics,\n".
        " Methods in Molecular Biology Volume 1231, 203-232, edited by A Mengoni, M Galardini and M Fondi.\n";
  print "\nThis software employs code, binaries and data from different authors, please cite them accordingly:\n";
  print " OrthoMCL v1.4 (www.orthomcl.org , PubMed:12952885)\n";
  print " INPARANOID v3.0 (inparanoid.sbc.su.se , PubMed:16873526)\n" if(!$DISSABLEPRND);
  print " COGtriangles v2.1 (sourceforge.net/projects/cogtriangles , PubMed=20439257)\n";
  print " NCBI Blast-2.2 (blast.ncbi.nlm.nih.gov , PubMed=9254694,20003500)\n";
  print " Bioperl v1.5.2 (www.bioperl.org , PubMed=12368254)\n";
  print " HMMER 3.1b2 (hmmer.org)\n";
  print " Pfam (http://pfam.xfam.org , PubMed=26673716)\n";
  print " Diamond v0.8.25 (https://github.com/bbuchfink/diamond, PubMed=25402007)\n";

  # check all binaries and data needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'d'}))
{
  $inputDIR = $opts{'d'};
  die "# EXIT : need a valid input directory\n" if(!-e $inputDIR);
  if(basename($inputDIR) =~ /([\+])/)
  {
    die "# EXIT : need a valid input directory name, offending char: '$1'\n";
  }
  if(substr($inputDIR,length($inputDIR)-1,1) eq '/'){ chop $inputDIR }
  $newDIR = $pwd.basename($inputDIR)."_homologues";
}
elsif(defined($opts{'i'}))
{
  $input_FASTA_file = $opts{'i'};
  die "# EXIT : need a valid input FASTA file\n" if(!-e $input_FASTA_file);
  my $root = basename($input_FASTA_file);
  if($root =~ /([\+])/)
  {
    die "# EXIT : need a valid input FASTA file, offending char: '$1'\n";
  }
  while($root =~ /\./){ $root = substr($root,0,-1) }
  $newDIR = $pwd.$root."_homologues";
}
else{ die "# EXIT : need a -d directory with .faa/.gbk files or a -i .faa file as input\n"; }

#if($input_FASTA_file){ $runmode = 'local'; }
if(defined($opts{'m'}))
{
  $runmode = $opts{'m'};
  if($runmode eq 'pbs'){ $runmode = 'cluster'; } # legacy
  elsif($runmode ne 'local' && $runmode ne 'cluster'){ $runmode = 'cluster'; }# default
}
else{ $runmode = 'local'; }

if($runmode eq 'local' && defined($opts{'n'}) && $opts{'n'} > 0)
{
  $n_of_cpus = $opts{'n'};
  $BLAST_NOCPU = $n_of_cpus;
}

if($runmode eq 'cluster' && !cluster_is_available())
{
  print "# EXIT : SGE cluster is not available, please appropriately set variable \$SGEPATH ($SGEPATH)\n";
  die "# EXIT : or choose runmode -m local\n";
}

if(defined($opts{'o'})){ $onlyblast = 1; }

if($opts{'r'}){ $reference_proteome_string = $opts{'r'}; }
else{ $reference_proteome_string = 0; }

check_installed_features(@FEATURES2CHECK);

if(defined($opts{'X'}))
{ 
  if(feature_is_installed('DIAMOND'))
  {
    $do_diamond = 1;
    $output_mask .= "dmd_";
    $pancore_mask .= "_dmd";
  }
  else{ warn_missing_soft('DIAMOND') }  
}

if(defined($opts{'f'}) && $opts{'f'} > 0 && $opts{'f'} <= 100)
{
  $filter_by_length = $opts{'f'};
  $output_mask .= "f$filter_by_length\_";
}
else{ $filter_by_length = 0; $output_mask .= "f0_"; }

if(defined($opts{'t'}) && $opts{'t'} >= 0)
{
  $min_cluster_size = $opts{'t'};
  $output_mask .= $min_cluster_size."taxa_";
}
else{ $min_cluster_size = 'all'; $output_mask .= "alltaxa_"; }

if($opts{'I'} && $inputDIR)
{
  $include_file = $opts{'I'};
  $output_mask .= basename($include_file)."_";
  $pancore_mask = "_".basename($include_file);
}

if($opts{'a'} && defined($opts{'d'}))
{
  if($do_diamond)
  {
    die "# EXIT : cannot run diamond (-X) with DNA features (-a)\n";
  }

  $do_features = $opts{'a'};
  $do_features =~ s/\s+//g;

  my $features = $do_features;
  $features =~ s/,//g;
  $output_mask .= "$features\_";
  $pancore_mask .= "_$features";
  rename_blast_outfiles('features');
}
else{ $do_features = 0 }

if(defined($opts{'A'})){ $do_ANIb_matrix = 1 }

if(defined($opts{'g'}) && defined($opts{'d'}) && !$do_features){ $do_intergenic = 1; }
else{ $do_intergenic = 0 }

if(defined($opts{'M'}))
{
  if(feature_is_installed('OMCL'))
  {
    $doMCL = 1;
    $output_mask .= "algOMCL_";
    $pancore_mask .= "_algOMCL";
  }
  else{ warn_missing_soft('OMCL') }
}
elsif(defined($opts{'P'}) && !$DISSABLEPRND)
{
  if(feature_is_installed('PRND'))
  {
    $doPARANOID = 1;
    $output_mask .= "algPRND_";
    $pancore_mask .= "_algPRND";
  }
  else{ warn_missing_soft('PRND') }
}
elsif(defined($opts{'G'})) # && $inputDIR ne '')
{
  if(feature_is_installed('COGS'))
  {
    $doCOG = 1;
    if(defined($opts{'x'}))
    {
      $COGmulticluster = 1;
      $output_mask .= "algCOGx_";
      $pancore_mask .= "_algCOGx";
    }
    else
    {
      $COGmulticluster = 0;
      $output_mask .= "algCOG_";
      $pancore_mask .= "_algCOG";
    }
  }
  else{ warn_missing_soft('COGS') }
}
else
{
  if(defined($opts{'b'}))
  {
    if($do_diamond)
    {
      die "# EXIT : cannot run diamond (-X) with -b option\n";
    }
    $do_minimal_BDBHs = 1;
    $output_mask .= "algBDBHmin_";
    $pancore_mask .= "_algBDBHmin";
  }
  else
  {
    $output_mask .= "algBDBH_";
    $pancore_mask .= "_algBDBH";
  }

  #if(defined($opts{'D'}) && !$do_features){}

  if($min_cluster_size eq '0')
  {
    die "\n# WARNING: use of the default BDBH algorithm with option -t 0 is not supported ".
      "(please check the manual)\n\n";
  }
  elsif($do_ANIb_matrix)
  {
    die "\n# WARNING: use of the default BDBH algorithm with option -A is not supported ".
      "(please check the manual)\n\n";
  }
}

if(defined($opts{'D'}) && !$do_features)
{
  if(feature_is_installed('PFAM'))
  {
    $do_PFAM = 1;
    $output_mask .= "Pfam\_";
    $pancore_mask .= "_Pfam";

    # in local mode increase batch size before calling _split_hmmscan.pl
    if($runmode eq 'local'){ $MAXPFAMSEQS = 5000 } # should fit most bacterial proteomes
  }
  else{ warn_missing_soft('PFAM') }
}

if(defined($opts{'s'}) && eval{ require DB_File } ){ $saveRAM = 1; }
else{ $saveRAM = 0; }

if(defined($opts{'c'}) && !$do_minimal_BDBHs)
{
  $do_genome_composition = 1;
  
  if(defined($opts{'z'}) && ($doMCL || $doCOG))
  {
    $do_soft = 1;
  }
  
  if($opts{'R'})
  {
    $random_number_generator_seed = $opts{'R'};
  }
  elsif($do_features)
  {
    die "\n# WARNING: use of -c and -a options usually requires setting -R ".
      "(please check the manual)\n\n";
  }
  
  srand($random_number_generator_seed);
}
else{ $do_genome_composition = 0; }

if(defined($opts{'e'}))
{
  $exclude_inparalogues = 1;
  $output_mask .= "e1_";
}
else{ $output_mask .= "e0_"; }

if(defined($opts{'E'}))
{
  $evalue_cutoff = $opts{'E'};
  if($evalue_cutoff > $MAXEVALUEBLASTSEARCH){ $evalue_cutoff = $MAXEVALUEBLASTSEARCH }
  $output_mask .= "E$evalue_cutoff\_"; $pancore_mask .= "_E$evalue_cutoff";
}
if(defined($opts{'C'}))
{
  $pmatch_cutoff = $opts{'C'}; # BDBH|OMCL
  if($pmatch_cutoff < 1){ $pmatch_cutoff = 1 }
  elsif($pmatch_cutoff > 100){ $pmatch_cutoff = 100 }
  $segment_cover_cutoff = $pmatch_cutoff; # COGS|PRND
  $output_mask .= "C$pmatch_cutoff\_"; $pancore_mask .= "_C$pmatch_cutoff";
}
if(defined($opts{'S'}))
{
  $pi_cutoff = $opts{'S'};
  if($pi_cutoff < 1){ $pi_cutoff = 1 }
  elsif($pi_cutoff > 100){ $pi_cutoff = 100 }
  $output_mask .= "S$pi_cutoff\_"; $pancore_mask .= "_S$pi_cutoff";
}
if(defined($opts{'F'}))
{
  $MCLinflation = $opts{'F'};
  if($MCLinflation < 1){ $MCLinflation = 1 }
  elsif($MCLinflation > 5){ $MCLinflation = 5 }
  $output_mask .= "F$MCLinflation\_"; $pancore_mask .= "_F$MCLinflation";
}
if(defined($opts{'N'}))
{
  $neighbor_corr_cutoff = $opts{'N'};
  if($neighbor_corr_cutoff < 0){ $neighbor_corr_cutoff = 0 }
  elsif($neighbor_corr_cutoff > 1){ $neighbor_corr_cutoff = 1 }
  $output_mask .= "N$neighbor_corr_cutoff\_"; $pancore_mask .= "_N$neighbor_corr_cutoff";
}
if(defined($opts{'B'}) && !$DISSABLEPRND)
{
  $bitscore_cutoff = $opts{'B'};
  if($bitscore_cutoff < 0){ $bitscore_cutoff = 0 }
  $output_mask .= "B$bitscore_cutoff\_"; $pancore_mask .= "_B$bitscore_cutoff";
}

print "# $0 -i $input_FASTA_file -d $inputDIR -o $onlyblast -X $do_diamond -e $exclude_inparalogues -f $filter_by_length -r $reference_proteome_string ".
  "-t $min_cluster_size -c $do_genome_composition -z $do_soft -I $include_file -m $runmode -n $n_of_cpus -M $doMCL -G $doCOG -P $doPARANOID ".
  "-C $pmatch_cutoff -S $pi_cutoff -E $evalue_cutoff -F $MCLinflation -N $neighbor_corr_cutoff -B $bitscore_cutoff -b $do_minimal_BDBHs ".
  "-s $saveRAM -D $do_PFAM -g $do_intergenic -a '$do_features' -x $COGmulticluster -R $random_number_generator_seed -A $do_ANIb_matrix\n\n";
  
if($runmode eq 'cluster')
{
  print "# computer cluster settings: $SGEPATH , $QUEUESETTINGS , $QUEUEWAIT , $WAITTIME\n\n";
}

###############################################################

## 0) declare most important vars
my ($infile,$new_infile,$p2oinfile,$seq,$comma_input_files,@newfiles,%ressize,%intergenic_FNA_files);
my ($label,%orthologues,$gene,$orth,$para,%inparalogues,%paralogues,$FASTAresultsDIR,$order,%orth_taxa,$minlog);
my ($n_of_similar_length_orthologues,$clusterfile,$dna_clusterfile,$previous_files,$current_files,$inpara);
my ($min_proteome_size,$reference_proteome,$smallest_proteome,$proteome_size,%seq_length,$cluster,$annot);
my ($smallest_proteome_name,$reference_proteome_name,%psize,$taxon,$n_of_clusters,$n_of_taxa,$n_of_residues);
my ($pname,$n_of_sequences,$refOK,$genbankOK,$cluster_size,$dna_cluster_size,$dna_cluster,$prot_cluster);
my (%orth_registry,%LSE_registry,$LSE_reference,$LSE_t,$redo_inp,$redo_orth,%idclusters,%names,$generef);
my ($reparse_all,$total_genbankOK,$n_of_similar_length_paralogues,$pfam_annot,$dnaOK,$n_of_taxa_cluster) = (0,0);
my ($BDBHdone,$PARANOIDdone,$orthoMCLdone,$COGdone,$n_of_parsed_lines,$n_of_pfam_parsed_lines) = (0,0,0,0,0,0);
my ($diff_BDBH_params,$diff_INP_params,$diff_HOM_params,$diff_OMCL_params,$lockcapableFS) = (0,0,0,0,0);
my ($total_clustersOK,$clgene,$clorth,%ANIb_matrix,%GIclusters,$clustersOK,%cluster_ids) = (0);

constructDirectory($newDIR);

# 0.1) try to make sure there is only 1 instance writing to $newDIR
# http://www.perlmonks.org/?node_id=590619, flock might fail with NFS filehandles
# http://stackoverflow.com/questions/4502721/in-perl-how-can-i-check-if-a-file-is-locked
# http://perldoc.perl.org/File/Temp.html
my ($fhtest,$testlockfilename) = tempfile( DIR => $newDIR );
if(flock($fhtest,LOCK_EX|LOCK_NB))
{
  $lockcapableFS = 1;
}
else
{
  print "# WARNING : cannot lock files in $newDIR ,\n".
    "# please ensure that no other instance of the program is running at this location\n\n";
}
unlink($testlockfilename);

open(my $fhlock,">$marfil_homology::lockfile") ||
  die "# EXIT : cannot create lockfile $marfil_homology::lockfile\n";
if($lockcapableFS)
{
  flock($fhlock, LOCK_EX|LOCK_NB) ||
    die "# EXIT : cannot run another instance of the program with same input data while previous is running\n";
}

# 0.2) prepare to tie hash/array to db if required : allows massives,slow data structs in disk
my (@sequence_data,@sequence_prot,@sequence_dna);
my $sequence_data_filename = $newDIR."/sequence_data.db";
my $sequence_prot_filename = $newDIR."/sequence_prot.db";
my $sequence_dna_filename = $newDIR."/sequence_dna.db";
my $pfam_data_filename = $newDIR."/pfam.db";# %pfam_hash is global, imported from marfil_homology
unlink($sequence_data_filename,$sequence_prot_filename,$sequence_dna_filename,$pfam_data_filename);
my $input_order_file = $newDIR."/input_order.txt";
my $diamond_file = $newDIR."/diamond.txt";
my $paranoidDIR = $newDIR."/paranoid_sql_C$pmatch_cutoff\_B$bitscore_cutoff\_O$segment_cover_cutoff/";

print "# version $VERSION\n";
print "# results_directory=$newDIR\n";
print "# parameters: MAXEVALUEBLASTSEARCH=$MAXEVALUEBLASTSEARCH MAXPFAMSEQS=$MAXPFAMSEQS ".
  "BATCHSIZE=$BATCHSIZE KEEPSCNDHSPS=$KEEPSCNDHSPS\n";
  
# make sure that DIAMOND jobs are not mixed with BLASTP jobs
if(-s $diamond_file)
{
  open(DMDFILE,$diamond_file) || die "# EXIT : cannot read $diamond_file\n";
  while(<DMDFILE>)
  {
    if(/^diamond job:(\d)/)
    {
      if($1 != $do_diamond)
      {
        $do_diamond = $1;
        print "# diamond job:$do_diamond (forced to match previous jobs)\n";
      }
      else{ print "# diamond job:$do_diamond\n" }  
    }    
  }
  close(DMDFILE);
}    
else
{ 
  if(-s $input_order_file)
  {
    # assume this was a BLASTP job to be sure
    if($do_diamond == 1)
    {
      $do_diamond = 0;
      print "# diamond job:$do_diamond (forced to match previous jobs)\n";
    }
    else{ print "# diamond job:$do_diamond\n" }
  }
  else{ print "# diamond job:$do_diamond\n" }
  
  open(DMDFILE,'>',$diamond_file) || die "# EXIT : cannot create $diamond_file\n";
  print DMDFILE "diamond job:$do_diamond\n";
  close(DMDFILE);
}

################################################################

## 1) read all input files, identify format and generate input FASTA files in temporary directory

print "\n# checking input files...\n";
$min_proteome_size = -1;
$reference_proteome = $refOK = $n_of_sequences = $n_of_taxa = $n_of_residues = 0;
$previous_files = $current_files = '';

# 1.0) instantiate BerkeleyDB if required
if($saveRAM)
{
  eval
  {
    import DB_File;
    tie(@sequence_data,'DB_File',$sequence_data_filename,1,0666,$DB_File::DB_RECNO)
      || die "# EXIT : cannot create $sequence_data_filename: $! (BerkeleyDB::Error)\n";
    tie(@sequence_prot,'DB_File',$sequence_prot_filename,1,0666,$DB_File::DB_RECNO)
      || die "# EXIT : cannot create $sequence_prot_filename: $! (BerkeleyDB::Error)\n";
    tie(@sequence_dna,'DB_File',$sequence_dna_filename,1,0666,$DB_File::DB_RECNO)
      || die "# EXIT : cannot create $sequence_dna_filename: $! (BerkeleyDB::Error)\n";

    if($do_PFAM && !$onlyblast)
    {
      tie(%pfam_hash,'DB_File',$pfam_data_filename,1,0666,$DB_File::DB_HASH)
        || die "# EXIT : cannot create $pfam_data_filename: $! (BerkeleyDB::Error)\n";
    }
  }
}

# 1.1) directory with input files
if($inputDIR)
{
  # 1.1.1) open and read directory, only master .fa/.gb files are considered
  opendir(DIR,$inputDIR) || die "# EXIT : cannot list $inputDIR\n";
  my @inputfiles = sort grep {/\.fa[a-z]*\S*$/i||/\.gb\S*$/||/_$/||/\.clusters$/} readdir(DIR);
  closedir(DIR);

  # 1.1.2) sort input files and put new files towards the end of @inputfiles: LILO
  if(-s $input_order_file)
  {
    my (@new_order_input_files,%previous_input_file,$n_of_new_infiles);

    open(ORDER,$input_order_file) || die "# EXIT : cannot read $input_order_file\n";
    while(<ORDER>)
    {
      chomp;
      ($order,$infile) = split(/\t/);
      if(!-s $inputDIR."/".$infile)
      { die "# EXIT : cannot find previous input file $infile, please re-run everything\n"; }

      $previous_input_file{$infile} = 1;
      $new_order_input_files[$order] = $infile;
    }
    close(ORDER);

    $n_of_new_infiles=0;
    foreach $infile (@inputfiles)
    {
      next if($previous_input_file{$infile});

      $new_order_input_files[++$order] = $infile;
      print "# order of new input file $infile = $order\n";
      $n_of_new_infiles++;
    }

    if($n_of_new_infiles){ print "# updating $input_order_file with $n_of_new_infiles new input files\n"; }
    open(ORDER,">$input_order_file") || die "# EXIT : cannot write $input_order_file\n";
    $order=0;
    foreach $infile (@new_order_input_files)
    {
      print ORDER "$order\t$infile\n";
      $order++;
    }
    close(ORDER);

    @inputfiles = @new_order_input_files;
  }
  else
  {
    $order=0;
    open(ORDER,">$input_order_file") || die "# EXIT : cannot write $input_order_file\n";
    foreach $infile (@inputfiles)
    {
      print ORDER "$order\t$infile\n";
      $order++;
    }
    close(ORDER);
  }

  # 1.1.3) check if all input files have GenBank format (-a/-g options)
  if($do_intergenic || $do_features)
  {
    foreach $infile (@inputfiles)
    {
      if($infile !~ /\.gb\S*$/i)
      {
        die "\n# EXIT : cannot analyze feature/intergenic sequences when input file is not in GenBank format\n\n";
      }
    }
  }

  # 1.1.4) open CSV protein to organism file, id to name (required by COGtriangles)
  open(COGCSV,">$p2ofilename") || die "# EXIT : cannot write $p2ofilename\n";

  # 1.1.5) iteratively parse input files
  foreach $infile (@inputfiles)
  {
    my ($file_ressize,$intergenic_size,$fasta_ref,$dna_fasta_ref) = (0,0);
    ($genbankOK,$dnaOK,$clustersOK,$proteome_size) = (0,0,0,-1);
    ++$n_of_taxa;

    # 1.1.5.1) genbank files
    if($infile =~ /\.gb\S*$/i)
    {
      my $prot_tmp_file = $newDIR ."/" .basename($infile) . ".faa";
      my $dna_tmp_file =  $newDIR ."/" .basename($infile) . ".fna";
      my $feature_dna_tmp_file = $newDIR ."/" .basename($infile) . "_features.fna";
      my $interg_tmp_file = $newDIR ."/" .basename($infile) .
        "_intergenic$MININTERGENESIZE\_$MAXINTERGENESIZE\_$INTERGENEFLANKORF.fna";

      if($do_features)
      {
        if(-s $feature_dna_tmp_file ||
          extract_features_from_genbank($inputDIR."/".$infile,
            $feature_dna_tmp_file,$do_features) > 0)
        {
          $fasta_ref = read_FASTA_file_array($feature_dna_tmp_file);
          if(!@$fasta_ref)
          {
            die "# ERROR: could not extract feature sequences from file ".
              "$inputDIR/$infile\n";
          }
        }
        else{ die "# ERROR: could not extract feature sequences from file $inputDIR/$infile\n"; }
      }
      else
      {
        if(-s $prot_tmp_file ||
          extract_CDSs_from_genbank($inputDIR."/".$infile,$prot_tmp_file,$dna_tmp_file) > 0)
        {
          $fasta_ref = read_FASTA_file_array($prot_tmp_file);
          if(-s $dna_tmp_file && @$fasta_ref)
          {
            $dna_fasta_ref = read_FASTA_file_array($dna_tmp_file);
            $genbankOK = 1;
            $total_genbankOK++;
          }
          else
          {
            print "# WARNING: could not extract nucleotide sequences from file $inputDIR/$infile\n";
          }
        }
        elsif(extract_genes_from_genbank($inputDIR."/".$infile,$prot_tmp_file,$dna_tmp_file) > 0)
        {
          print "# WARNING: can only extract genes (not CDSs) from file $inputDIR/$infile\n";

          # both dna & prot sequences are always generated with extract_genes
          $fasta_ref = read_FASTA_file_array($prot_tmp_file);
          $dna_fasta_ref = read_FASTA_file_array($dna_tmp_file);
          $genbankOK = 1;
          $total_genbankOK++;
        }
        else{ die "# ERROR: could not extract sequences from file $inputDIR/$infile\n"; }

        if($do_intergenic)
        {
          if(!-s $interg_tmp_file)
          {
            $intergenic_size =
              extract_intergenic_from_genbank($inputDIR."/".$infile,
              $interg_tmp_file,$MININTERGENESIZE,
              $MAXINTERGENESIZE,$INTERGENEFLANKORF);
          }
          else
          {
            my $intergfna_ref = read_FASTA_file_array($interg_tmp_file);
            $intergenic_size = scalar(@$intergfna_ref);
          }

          if(!$intergenic_size)
          {
            print "# WARNING: could not extract intergenic sequences from file $inputDIR/$infile\n";
          }
          else{ $intergenic_FNA_files{$infile} = $interg_tmp_file }
        }
      }
    }
    elsif(($infile =~ /\.clusters$/ || $infile =~ /_$/) && -d $inputDIR.'/'.$infile) # 1.1.5.2) subdirectory with cluster files in FASTA format
    {
      if($total_clustersOK++)
      {
        die "\n# WARNING: can only take one folder of pre-clustered sequences\n\n";
      }
      
      
      opendir(CLDIR,$inputDIR.'/'.$infile) || die "# EXIT : cannot list $inputDIR/$infile\n";
      my @clusterfiles = sort grep {/\.fa[a-z]*$/i} readdir(CLDIR);
      closedir(CLDIR);

      my ($cluster_id,$cldnaOK) = (0);
      my $clseq = $n_of_sequences + 1; # make sure sequences in clusters are correctly numbered
      my (@clusters_fasta,@clusters_dna_fasta);

      foreach my $clusterfile (@clusterfiles)
      {
        $cldnaOK = 0;
        my $protfile = $inputDIR.'/'.$infile.'/'.$clusterfile; 
        $fasta_ref = read_FASTA_file_array($protfile);
        
        # check whether an identical .fna file exists (expects sequences in same order)
        my $dnafile = $protfile;
        $dnafile =~ s/\.fa+.*?$/\.fna/;
        if($dnafile ne $protfile && -s $dnafile )
        {
          $dna_fasta_ref = read_FASTA_file_array($dnafile);
          if($#{$fasta_ref} == $#{$dna_fasta_ref}){ $cldnaOK = 1 }
          else
          {
            print "WARNING: cannot use nucleotide sequences in file $dnafile\n".
              "# as they do not match those in file $protfile\n";
          }
        }
        
        # arbitrarily choose first cluster sequence as representative
        $cluster_id = $clseq;
        
        # add sequentially members of pre-processed cluster
        foreach $seq ( 0 .. $#{$fasta_ref} )
        {
          # remember members of this cluster
          push(@{$cluster_ids{$cluster_id}},$clseq) if($clseq != $cluster_id); 
        
          $clusters_fasta[$clseq] = $fasta_ref->[$seq]; #print "$clseq\n";
          if($cldnaOK)
          {
            $clusters_dna_fasta[$clseq] = $dna_fasta_ref->[$seq]; 
          }  
          $clseq++;
        }
      }
      
      $fasta_ref = \@clusters_fasta;
      $dna_fasta_ref = \@clusters_dna_fasta;
      
      if($#{$fasta_ref} == $#{$dna_fasta_ref}){ $dnaOK = 1 }
      $clustersOK = scalar(@clusterfiles);
      $total_clustersOK++;
      
      if($do_genome_composition)
      {
        die "\n# WARNING: -c genome composition analyses cannot be done with pre-clustered input sequences ".
          "(please check the manual)\n\n";
      }
    }
    else # 1.1.5.3) fasta files
    {
      my $protfile = $inputDIR."/".$infile;
      $fasta_ref = read_FASTA_file_array($protfile);

      # check whether an identical .fna file exists (expects sequences in same order)
      my $dnafile = $protfile;
      $dnafile =~ s/\.fa+.*?$/\.fna/;
      if($dnafile ne $protfile && -s $dnafile )
      {
        $dna_fasta_ref = read_FASTA_file_array($dnafile);
        if($#{$fasta_ref} == $#{$dna_fasta_ref}){ $dnaOK = 1 }
        else
        {
          print "WARNING: cannot use nucleotide sequences in file $dnafile\n".
            "# as they do not match those in file $protfile\n";
        }
      }
    }

    if(!defined($fasta_ref)){ die "# ERROR: could not extract sequences from $inputDIR/$infile\n"; }
    $proteome_size = scalar(@$fasta_ref);
    $psize{$infile} = $proteome_size;

    if($do_intergenic){ print "# $infile $proteome_size (intergenes=$intergenic_size)\n"; }
    elsif($clustersOK)
    { 
      if($dnaOK){ print "# $infile $proteome_size ($clustersOK clusters, found twin .fna file)\n"; }
      else{ print "# $infile $proteome_size ($clustersOK clusters)\n"; }
    }
    elsif($dnaOK){ print "# $infile $proteome_size (found twin .fna file)\n"; }
    else{ print "# $infile $proteome_size\n"; }

    $new_infile = basename($infile) . ".$FASTAEXTENSION";
    $p2oinfile  = basename($infile);
    if($do_features)
    {
      $new_infile = basename($infile) . "_features.$FASTAEXTENSION";
      $p2oinfile  = basename($infile) . '_features';
    }

    open(INPUT,">$newDIR/$new_infile") || die "# EXIT : cannot create $newDIR/$new_infile\n";
    foreach $seq ( 0 .. $#{$fasta_ref} )
    {
      if(!$fasta_ref->[$seq][SEQ] || $fasta_ref->[$seq][SEQ] eq '')
      {
        chomp $fasta_ref->[$seq][NAME];
        print "# skipped void sequence ($fasta_ref->[$seq][NAME])\n" if($fasta_ref->[$seq][NAME]);
        next;
      }

      # skip DNA sequences in supposed-to-be amino acid files
      if(!$do_features && $fasta_ref->[$seq][SEQ] =~ /^[ACGTUXN\-\s]+$/i)
      {
        print "# skipped nucleotide sequence ($fasta_ref->[$seq][NAME])\n";
        next;
      }

      # @sequence_... global arrays take large amounts of RAM
      # first sequence is [1]
      ++$n_of_sequences;
      $file_ressize += length($fasta_ref->[$seq][SEQ]);
      if(!$onlyblast)
      {
        $sequence_data[$n_of_sequences] = $fasta_ref->[$seq][NAME];
        $sequence_prot[$n_of_sequences] = $fasta_ref->[$seq][SEQ]; # with -a actually stores dna
        if($fasta_ref->[$seq][SEQ] eq '')
        {
          printf("# cannot parse protein sequence of %s... (%s)\n",
            substr($fasta_ref->[$seq][NAME],0,15),$infile);
        }

        if($genbankOK || $dnaOK)
        {
          $sequence_dna[$n_of_sequences] = $dna_fasta_ref->[$seq][SEQ];
          if($dna_fasta_ref->[$seq][SEQ] eq '')
          {
            printf("# cannot parse DNA sequence of %s... (%s)\n",
              substr($fasta_ref->[$seq][NAME],0,15),$infile);
          }
        }
      }
      
      # do not print redundant pre-clustered sequences, take only 1/cluster
      next if($clustersOK && !$cluster_ids{$seq});
      
      print INPUT ">$n_of_sequences\n$fasta_ref->[$seq][SEQ]\n";
      $pname = $fasta_ref->[$seq][NAME];$pname =~ s/#/hash/g; # COG does not like hashes in names
      print COGCSV "$n_of_sequences,$p2oinfile,$pname\n";
    }
    close(INPUT);

    $comma_input_files .= "$new_infile,";
    push(@newfiles,$new_infile);

    $n_of_residues += $file_ressize;
    $ressize{$infile} = $ressize{$new_infile} = $file_ressize;
  }

  close(COGCSV);

  if($min_cluster_size eq 'all'){ $min_cluster_size = $n_of_taxa; } # size of clusters in section 4
}
else # 1.2) input single FASTA file with [taxon names] in headers, order not checked
{
  print "# $input_FASTA_file\n";

  # 1.2.0) open CSV protein to organism file required by COGtriangles
  open(COGCSV,">$p2ofilename") || die "# EXIT : cannot write $p2ofilename\n";

  my $fasta_ref = read_FASTA_file_array($input_FASTA_file);

  # allow taxon/genome name to be a single word (1)
  my %taxa_in_headers = find_taxa_FASTA_array_headers($fasta_ref,1);

  # 1.2.1) print sequences of the same taxon in separate files
  foreach my $taxon (sort (keys(%taxa_in_headers)))
  {
    my $taxon_ressize = 0;
    $proteome_size = 0;

    ++$n_of_taxa;

    $taxon =~ s/\[|\]//g;
    $new_infile = $taxon;
    $taxon =~ s/\(/\\(/g;$taxon =~ s/\)/\\)/g; # avoid problems with headers containing ()s
    $new_infile =~ s/\s+|\/|\.|\||\/|\*|\;|\,|"|'|\(|\)/_/g;
    $new_infile .= "_$n_of_taxa.faa"; #print "$new_infile\n";

    open(INPUT,">$newDIR/$new_infile.$FASTAEXTENSION")
      || die "# EXIT : cannot create $newDIR/$new_infile.$FASTAEXTENSION\n";
    foreach $seq ( 0 .. $#{$fasta_ref} )
    {
      if($fasta_ref->[$seq][NAME] =~ /$taxon/)
      {
        # skip DNA-only sequences
        if($fasta_ref->[$seq][SEQ] =~ /^[ACGTUXN\-\s]+$/i)
        {
          print "# skipped nucleotide sequence ($fasta_ref->[$seq][NAME])\n";
          next;
        }

        # @sequence_... global arrays take large amounts of RAM
        # first sequence is [1]
        $n_of_sequences++;
        $taxon_ressize += length($fasta_ref->[$seq][SEQ]);
        if(!$onlyblast)
        {
          $sequence_data[$n_of_sequences] = $fasta_ref->[$seq][NAME];
          $sequence_prot[$n_of_sequences] = $fasta_ref->[$seq][SEQ];
        }

        print INPUT ">$n_of_sequences\n$fasta_ref->[$seq][SEQ]\n";
        $pname = $fasta_ref->[$seq][SEQ];$pname =~ s/#/hash/g; # COG does not like hashes in names
        print COGCSV "$n_of_sequences,$new_infile.$FASTAEXTENSION,$pname\n";
        $proteome_size++;
      }
    }
    close(INPUT);
    $comma_input_files .= "$new_infile.$FASTAEXTENSION,"; # undesired genomes should not be added to this string
    push(@newfiles,"$new_infile.$FASTAEXTENSION");
    $psize{$new_infile} = $proteome_size;

    $n_of_residues += $taxon_ressize;
    $ressize{$new_infile} = $taxon_ressize;
  }

  close(COGCSV);

  if($min_cluster_size eq 'all'){ $min_cluster_size = $n_of_taxa; } # size of clusters in section 4
}

# 1.2.5) read taxa labels and concatenate FASTA files ($all_fa_file is put in $newDIR/tmp, global @taxa is filled in here)
chop($comma_input_files);
if($inputDIR)
{ 
  %seq_length = %{ constructAllFasta($comma_input_files,$n_of_sequences,'',$filter_by_length) }; 
}
else{ %seq_length = %{ constructAllFasta($comma_input_files,$n_of_sequences,$all_fa_file,$filter_by_length) } }

# 1.4) correct list of included files, must be done after reading them all for correct numbering
if($include_file)
{
  my ($included,$includedfull,$n_of_matched_included,@Itaxa);

  open(INCL,$include_file) || die "# EXIT : cannot read $include_file\n";
  while(<INCL>)
  {
    next if(/^#/ || /^$/);
    $included_input_files{(split)[0]} = $.;
  }
  close(INCL);

  print "# included input files (".scalar(keys(%included_input_files))."):\n";

  $refOK = $n_of_sequences = $n_of_residues = $n_of_matched_included = 0;
  TAXON: foreach $included (sort {$included_input_files{$a}<=>$included_input_files{$b}}
    keys(%included_input_files))
  {
    $includedfull = $included;
    if($do_features){ $includedfull .= '_features'; }
    foreach $taxon (@taxa)
    {
      if($includedfull =~ /$taxon/)
      {
        push(@Itaxa,$taxon);
        $n_of_sequences += $psize{$included};
        $n_of_residues  += $ressize{$included};
        print ": $taxon $included $psize{$included}\n";
        $n_of_matched_included++;
        next TAXON;
      }
    }
  }
  print "\n";

  if($n_of_matched_included < 3)
  {
    die "# EXIT : failed to match taxa included in $include_file ($n_of_matched_included), ".
      "please make sure their names match those of input files\n";
  }

  # update @taxa, $min_cluster_size & $total_genbankOK
  if($total_genbankOK == $n_of_taxa){ $total_genbankOK = scalar(@Itaxa) }
  else{} # perhaps some warning should go in here?
  @taxa = @Itaxa;
  $n_of_taxa = scalar(@Itaxa);
  if($min_cluster_size eq 'all'){ $min_cluster_size = $n_of_taxa; }
  elsif($min_cluster_size > $n_of_taxa){ $min_cluster_size = $n_of_taxa; }
}

# 1.4.1) set minimum bit-score for neighbor correlation (similarity) calculations (Song,Joseph,Davis,Durand(2008)PlOS PMID:18475320)
$MIN_BITSCORE_SIM_DEFAULT = sprintf("%1.3f", 0.50 * ( log(($n_of_residues ** 2) / ($n_of_sequences * $MAXEVALUEBLASTSEARCH) ) / log(2)) );

printf("# taxa considered = %d sequences = %d residues = %d MIN_BITSCORE_SIM = %1.1f\n\n",$n_of_taxa,$n_of_sequences,$n_of_residues,$MIN_BITSCORE_SIM_DEFAULT);

if($n_of_taxa<2){ die "# EXIT: need at least two taxa to make clusters\n" }
elsif($doCOG && $n_of_taxa<3){ die "# EXIT: need at least three taxa to run COGtriangles algorithm\n" }

#my $RAM = estimate_RAM($n_of_sequences);
#printf("# estimated memory requirements: %s\n\n",$RAM) if(!$onlyblast && $RAM);

# 1.5) set reference proteome index and mask (by default, select proteome with smallest number of sequences)
# it is necessary for BDBH in particular and for all algorithms when handling intergenic regions
for($taxon=0;$taxon<scalar(@taxa);$taxon++)
{
  $proteome_size = $psize{$taxa[$taxon]};

  # update minimal proteome size
  if($min_proteome_size == -1 || $proteome_size < $min_proteome_size)
  {
    $min_proteome_size = $proteome_size;
    $smallest_proteome_name = $taxa[$taxon];
    $smallest_proteome = $taxon;
  }

  # check user-defined reference proteome
  if($reference_proteome_string ne '0' &&
    !$refOK && $taxa[$taxon] =~ /$reference_proteome_string/)
  {
    $reference_proteome = $taxon;
    $reference_proteome_name = $taxa[$taxon];
    $refOK=1;
  }
}

if(!$refOK && $reference_proteome_string ne '0')
{
  print "# WARNING: cannot find reference proteome ($reference_proteome_string), taking default\n";
}

if($reference_proteome_string eq '' || !$refOK)
{
  $reference_proteome_name = $smallest_proteome_name;
  $reference_proteome = $smallest_proteome;
}
$reference_proteome_name = (split(/\./,$reference_proteome_name))[0];
$reference_proteome_name =~ s/[\s+|_]//g;
$output_mask = $reference_proteome_name."\_" . $output_mask;

print "# mask=$output_mask ($pancore_mask)\n" if(!$onlyblast);

# 1.6) check previously processed input files to decide whether blast parsing
# is needed again and update $selected_genomes_file
if($do_minimal_BDBHs){ $minlog = "minimal_BDBHs$reference_proteome" } # reference is crucial for -b runs
else{ $minlog = 'minimal_BDBHsNO' }
$current_files = join('',sort(@taxa,$minlog));
if(-s $selected_genomes_file)
{
  $previous_files = get_string_with_previous_genomes($selected_genomes_file);
  if($current_files eq $previous_files)
  {
    print "\n# skipped genome parsing (".short_path($selected_genomes_file,$pwd).")\n\n";
  } #else{ print "$current_files ne $previous_files\n" }
}

open(SEL,">$selected_genomes_file") || die "# cannot create $selected_genomes_file\n";
foreach $taxon (@taxa){ print SEL "$taxon\n"; }
print SEL "$minlog\n"; #
close(SEL);

# 1.7) run Pfam searches if requested (recommended only in cluster mode, otherwise it might get very slow)
if($do_PFAM)
{
  my ($pfamout,$hmmerout,$total_seqs_group,$group,$group_fasta_file,$group_clusterfile,$command,@group_pfamfiles);
  my ($group_pfamfile,$totalseqs,%cluster_PIDs,%pfamgroups);

  print "\n# submitting Pfam HMMER jobs ... \n"; #(MAXPFAMSEQS=$MAXPFAMSEQS) ...\n";
  foreach $new_infile (@newfiles)
  {
    next if(-s $newDIR ."/_". $new_infile .".pfam"); # comprueba que no exista ya
    $infile = (split(/\.$FASTAEXTENSION/,$new_infile))[0];
    next if($include_file && !$included_input_files{$infile});

    # split Pfam searches in groups of $MAXPFAMSEQS searches
    my $fasta_ref = read_FASTA_file_array($newDIR ."/". $new_infile); # minimize RAM usage, file by file
    $totalseqs = scalar(@$fasta_ref);
    $total_seqs_group = $group = 0;
    $group_fasta_file = $newDIR . "/_". $new_infile . $group;
    $group_pfamfile   = $newDIR . "/_". $new_infile . $group . '.pfam';
    $group_clusterfile    = $newDIR . "/_". $new_infile . $group . '.queue';
    push(@group_pfamfiles,$group_fasta_file,$group_pfamfile,$group_clusterfile);

    open(FGR,">$group_fasta_file") || die "# EXIT : cannot create fasta group file $group_fasta_file : $!\n";

    foreach $seq (0 .. $#{$fasta_ref})
    {
      $total_seqs_group++;
      print FGR ">$fasta_ref->[$seq][NAME]\n$fasta_ref->[$seq][SEQ]\n";
      if($seq+1 == $totalseqs || ($MAXPFAMSEQS && $total_seqs_group == $MAXPFAMSEQS))
      {
        close(FGR);

        if(!-s $group_pfamfile) # scan Pfam only if necessary
        {
          $command = format_HMMPFAM_command()." $group_fasta_file > $group_pfamfile";

          if($runmode eq 'cluster')
          {
            submit_cluster_job($new_infile.$group,$command,$group_clusterfile,$newDIR,
              $SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
          }
          else # 'local' runmode
          {
            $command = format_SPLITHMMPFAM_command()."$BATCHSIZE $command ";
            system("$command");
            if($? != 0)
            {
              die "# EXIT: failed while running localPfam search ($command)\n";
            }
          }
        }
        $pfamgroups{$new_infile}{$group} = $group_pfamfile;

        if($seq < $totalseqs) # take care of remaining sequences/groups
        {
          $total_seqs_group = 0;
          $group++;
          $group_fasta_file = $newDIR . "/_". $new_infile . $group;
          $group_pfamfile   = $newDIR . "/_". $new_infile . $group . '.pfam';
          $group_clusterfile    = $newDIR . "/_". $new_infile . $group . '.queue';
          push(@group_pfamfiles,$group_fasta_file,$group_pfamfile,$group_clusterfile);
          open(FGR,">$group_fasta_file") || die "# $0 : cannot create fasta group file $group_fasta_file : $!\n";
        }
      }
    }
  }

  # wait until PBD Pfam jobs are done
  if($runmode eq 'cluster')
  {
    check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
  }
  print "# done\n\n";

  # concatenate partial/group Pfam results
  my (@pfamfiles);
  foreach $new_infile (@newfiles)
  {
    my @group_pfamfiles;
    $pfamout = $newDIR ."/_". $new_infile . ".pfam";
    push(@pfamfiles,$pfamout);
    next if(-s $pfamout);

    print "# concatenating Pfam files ($pfamout)...\n";
    foreach $group (sort {$a<=>$b} keys(%{$pfamgroups{$new_infile}}))
    {
      $group_pfamfile = $pfamgroups{$new_infile}{$group};
      if(!-e $group_pfamfile)
      {
        die "# EXIT, $group_pfamfile does not exist, Pfam search might have failed or ".
          "hard drive is still writing it (please re-run)\n";
      }
      system("cat $group_pfamfile >> $pfamout");
    }
    print "# done\n\n";
  }

  # remove group files if concatenation went OK, minimize disk usage
  unlink(@group_pfamfiles);

  # parse Pfam results if required (pfam_file is global)
  if(!-s $pfam_file || $current_files ne $previous_files) # || $include_file)
  {
    print "# parsing Pfam domain assignments (generating $pfam_file) ...\n\n";
    $n_of_pfam_parsed_lines = pfam_parse( $pfam_file , @pfamfiles );
  }
  else{ print "# re-using parsed Pfam assignment file ".short_path($pfam_file,$pwd)." ...\n\n"; }
}

################################################################

# 2) run BLAST/DIAMOND searches in created directory and parse output
if(!-s $bpo_file || $current_files ne $previous_files || ($doCOG && !-s $cogblasthits)) # || $include_file
{
  if($do_diamond)
  {
    if(!feature_is_installed('DIAMOND')){ warn_missing_soft('diamond blastp/diamond makedb') } 
  }
  else
  {
    if(!feature_is_installed('BLAST')){ warn_missing_soft('blastp/makeblastdb') } 
  }  

  my ($blastout,$clusteroutfile,%cluster_PIDs,@tmp_blast_output_files);
  my ($blastDBfile,$command,@to_be_deleted);

  # remove previous results to avoid confusion if required (-I,-a,new files)
  if($do_features || $current_files ne $previous_files){ unlink($blast_file,$bpo_file,$cogblasthits); } # $include_file ||

  if($do_diamond) # format own and rest FASTA files for DIAMOND
  {
    foreach $new_infile (@newfiles)
    {
      # format sequences of each individual accession (to put inparalogues up)
      if(!-s $newDIR ."/". $new_infile .".dmnd")
      {
        $infile = (split(/\.$FASTAEXTENSION/,$new_infile))[0];
        next if($include_file && !$included_input_files{$infile});
        executeDIAMONDMAKEDB($newDIR ."/". $new_infile);
        if(!-s $newDIR ."/". $new_infile .".dmnd")
        {
          die "# EXIT: cannot format DIAMOND protein sequence base $newDIR/$new_infile\n"
        }
      }
      
      # now concat the other accessions
      $blastDBfile = $newDIR ."/". $new_infile .".others";
      if(!-s $blastDBfile.".dmnd" || $current_files ne $previous_files)
      {
        open(ALLOTHERS,">",$blastDBfile) || die "# EXIT: cannot create $blastDBfile\n";
        foreach my $other_infile (@newfiles)
        {
          next if($other_infile eq $new_infile);
          $infile = (split(/\.$FASTAEXTENSION/,$other_infile))[0];
          next if($include_file && !$included_input_files{$infile}); 
          
          open(OTHER,"$newDIR/$other_infile") ||
            die "# EXIT: cannot read $newDIR/$other_infile\n";
          while(<OTHER>){ print ALLOTHERS }
          close(OTHER);  
        }  
        close(ALLOTHERS);  
          
        executeDIAMONDMAKEDB($blastDBfile);
        if(!-s $blastDBfile .".dmnd")
        {
          die "# EXIT: cannot format DIAMOND protein sequence base $blastDBfile\n"
        }
        else{ unlink($blastDBfile) } # we don't need this anymore
      }
    }
  }  
  else # format FASTA files for BLAST
  {
    foreach $new_infile (@newfiles)
    {
      if(!$do_features)
      {
        next if(-s $newDIR ."/". $new_infile .".psq");
        $infile = (split(/\.$FASTAEXTENSION/,$new_infile))[0];
        next if($include_file && !$included_input_files{$infile});
        
        executeFORMATDB($newDIR ."/".$new_infile);
        if(!-s $newDIR ."/". $new_infile .".psq")
        {
          die "# EXIT: cannot format BLAST protein sequence base $newDIR/$new_infile\n"
        }    
      }
      else
      {
        next if(-s $newDIR ."/". $new_infile .".nsq");
        $infile = (split(/_features\.$FASTAEXTENSION/,$new_infile))[0];
        next if($include_file && !$included_input_files{$infile});

        executeFORMATDB($newDIR ."/".$new_infile,1);
        if(!-s $newDIR ."/". $new_infile .".nsq")
        {
          die "# EXIT: cannot format BLAST nucleotide sequence base $newDIR/$new_infile\n"
        }
      }
    }
  }
   
  if($do_diamond) # launch DIAMOND jobs
  {
    print "\n# running DIAMOND searches ...\n";
    
    # self-scan each input accession
    foreach my $bfile1 (0 .. $#newfiles)
    {
      # self-scan each input accession
      $new_infile = $newfiles[$bfile1];
      $infile = (split(/\.$FASTAEXTENSION/,$new_infile))[0];
      next if($include_file && !$included_input_files{$infile});
      
      $blastDBfile = $newDIR ."/". $new_infile;
      $blastout   = $newDIR ."/_". $new_infile ."_". $new_infile.".blast";
      $clusteroutfile = $newDIR ."/_". $new_infile ."_". $new_infile.".queue";
      push(@to_be_deleted,$clusteroutfile);

      if(-s $blastout) # check previous runs
      {
        if(!-s $blast_file){ push(@tmp_blast_output_files,$blastout); }
        next;
      }
        
      $command = format_DIAMONDblastp_command("$newDIR/$new_infile",
        $blastout,$blastDBfile,$MAXEVALUEBLASTSEARCH,$psize{$infile});
            
      if($runmode eq 'cluster')
      {
        submit_cluster_job($new_infile,$command,$clusteroutfile,$newDIR,
          $SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
      }
      else # 'local' runmode -> ensure multi-core CPU use
      {
        # hack to use allocated number of threads in local mode
        $command =~ s/--threads 1/--threads $BLAST_NOCPU/;
        system("$command");
        if($? != 0)
        {
          die "# EXIT: failed while running local DIAMOND search ($command)\n";
        }
      }

      push(@tmp_blast_output_files,$blastout);    
    }
    
    # scan each accession against the others
    foreach my $bfile1 (0 .. $#newfiles)
    {
      $new_infile = $newfiles[$bfile1];
      $infile = (split(/\.$FASTAEXTENSION/,$new_infile))[0];
      next if($include_file && !$included_input_files{$infile});
      
      $blastDBfile = $newDIR ."/". $new_infile .".others";
      $blastout   = $newDIR ."/_". $new_infile ."_others.blast";
      $clusteroutfile = $newDIR ."/_". $new_infile ."_others.queue";
      push(@to_be_deleted,$clusteroutfile);

      if(-s $blastout && $current_files eq $previous_files) # check previous runs
      {
        if(!-s $blast_file){ push(@tmp_blast_output_files,$blastout); }
        next;
      }
        
      $command = format_DIAMONDblastp_command("$newDIR/$new_infile",
        $blastout,$blastDBfile,$MAXEVALUEBLASTSEARCH,$n_of_sequences);
            
      if($runmode eq 'cluster')
      {
        submit_cluster_job($new_infile,$command,$clusteroutfile,$newDIR,
          $SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
      }
      else # 'local' runmode -> ensure multi-core CPU use
      {
        $command =~ s/--threads 1/--threads $BLAST_NOCPU/;
        system("$command");
        if($? != 0)
        {
          die "# EXIT: failed while running local DIAMOND search ($command)\n";
        }
      }

      push(@tmp_blast_output_files,$blastout);     
    }         
  }
  else  # launch BLAST jobs
  {
    print "\n# running BLAST searches ...\n";
    foreach my $bfile1 (0 .. $#newfiles)
    {
      $new_infile = $newfiles[$bfile1];
      $infile = (split(/\.$FASTAEXTENSION/,$new_infile))[0];
      if($do_features){ $infile = (split(/_features\.$FASTAEXTENSION/,$new_infile))[0]; }
      next if($include_file && !$included_input_files{$infile});

      foreach my $bfile2 (0 .. $#newfiles)
      {
        $blastDBfile = $newfiles[$bfile2];
        $infile = (split(/\.$FASTAEXTENSION/,$blastDBfile))[0];
        if($do_features){ $infile = (split(/_features\.$FASTAEXTENSION/,$new_infile))[0]; }
        next if($include_file && !$included_input_files{$infile});

        $blastout   = $newDIR ."/_". $new_infile ."_". $blastDBfile.".blast";
        $clusteroutfile = $newDIR ."/_". $new_infile ."_". $blastDBfile.".queue";
        push(@to_be_deleted,$clusteroutfile);

        if($do_minimal_BDBHs) # skip non-minimal searches if required
        {
          next if( ($bfile1 != $reference_proteome && $bfile2 != $reference_proteome) &&
            ($bfile1 != $bfile2) );#print "$bfile1 $bfile2 $new_infile $blastDBfile\n";
        }
        if(-s $blastout) # check previous BLAST runs
        {
          if(!-s $blast_file){ push(@tmp_blast_output_files,$blastout); }
          next;
        }
        
        if($do_features)
        {
          $command = format_BLASTN_command("$newDIR/$new_infile",
            $blastout,"$newDIR/$blastDBfile",
            $MAXEVALUEBLASTSEARCH,$psize{$infile},$runmode ne 'cluster');
        }
        else
        {
          $command = format_BLASTP_command("$newDIR/$new_infile",
            $blastout,"$newDIR/$blastDBfile",
            $MAXEVALUEBLASTSEARCH,$psize{$infile},$runmode ne 'cluster');
        }   

        if($runmode eq 'cluster')
        {
          submit_cluster_job($new_infile,$command,$clusteroutfile,$newDIR,
            $SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
        }
        else # 'local' runmode -> ensure multi-core CPU use
        {
          $command = format_SPLITBLAST_command()."$BATCHSIZE $command > /dev/null"; 
          system("$command");
          if($? != 0)
          {
            die "# EXIT: failed while running local BLAST ($command)\n";
          }
        }

        push(@tmp_blast_output_files,$blastout);
      }
    }
  }

  # wait until cluster jobs are done
  if($runmode eq 'cluster')
  {
    check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
  }
  print "# done\n\n";

  # concatenate blast output files to $blast_file (global var)
  if(@tmp_blast_output_files)
  {
    print "# concatenating and sorting BLAST/DIAMOND results...\n";
    foreach $blastout (@tmp_blast_output_files)
    {
      if(!-e $blastout)
      {
        sleep($QUEUEWAIT); # give disk extra time
        if(!-e $blastout)
        {
          die "# EXIT, $blastout does not exist, BLAST/DIAMOND search might failed ".
            "or hard drive is still writing it (please re-run)\n";
        }
      }
    }
    
    # NxN BLAST jobs OR 2xN diamond jobs (for N input  files)
    sort_blast_results($blast_file,$KEEPSCNDHSPS,@tmp_blast_output_files);
      
    print "# done\n\n";
    unlink(@to_be_deleted); # remove .queue files
  }

  # parse search results and create $bpo_file
  if(!-s $bpo_file || $current_files ne $previous_files) # || $include_file)
  {
    $n_of_parsed_lines = blast_parse($blast_file,$bpo_file,\%seq_length,0,1);
    $reparse_all = 1;
  }

  if($doCOG && ($reparse_all || !-s $cogblasthits))
  {
    blast_parse_COG($blast_file,1);  # depends on $MAXEVALUEBLASTSEARCH
    $reparse_all = 1;
  }
}
else
{
  print "# skip BLAST searches and parsing\n\n".
    "# WARNING: please remove/rename results directory:\n# '$newDIR/'\n".
    "# if you change the sequences in your .gbk/.faa files or want to re-run\n";
}

undef(%seq_length) unless($filter_by_length); # free hash if not needed anymore

if($onlyblast)
{
  print "\n# terminating after BLAST (-o)\n";
  if($saveRAM)
  {
    untie(@sequence_data);untie(@sequence_prot);untie(@sequence_dna);
    unlink($sequence_data_filename,$sequence_prot_filename,$sequence_dna_filename);
  }

  my $end_time = new Benchmark();
  print "\n# runtime: ".timestr(timediff($end_time,$start_time),'all')."\n";
  print "# RAM use: ".calc_memory_footprint()."\n";
  exit(0);
}
else
{
  if($saveRAM)
  {
    print "\n# creating taxa indexes...\n";
    construct_taxa_indexes($bpo_file);
  }
  else
  {
    if(!$n_of_parsed_lines)
    {

      # count lines in $bpo_file in a portable manner, assuming \n
      #$n_of_parsed_lines = (split(/\s+/,`wc -l $bpo_file`))[0];
      open(BPO,$bpo_file) || die "# cannot open $bpo_file\n";
      $n_of_parsed_lines += tr/\n/\n/ while sysread(BPO,$_,2 ** 16);
      close(BPO);
    }
    printf("\n# creating indexes, this might take some time (lines=%1.2e) ...\n\n",$n_of_parsed_lines);
    if(!$n_of_parsed_lines)
    {
      die "# EXIT: parsed BLAST output ($bpo_file) seems to be empty, please remove '$newDIR/' and re-run\n";
    }
    construct_taxa_indexes($bpo_file); # needed for cluster nodes
    construct_indexes($bpo_file);
  }

  if($do_PFAM)
  {
    if(!$n_of_pfam_parsed_lines)
    {
      open(PFAM,$pfam_file) || die "# cannot open $pfam_file\n";
      $n_of_pfam_parsed_lines += tr/\n/\n/ while sysread(PFAM,$_,2 ** 16);
      close(PFAM); #$n_of_pfam_parsed_lines = (split(/\s+/,`wc -l $pfam_file`))[0];
    }

    printf("\n# creating Pfam indexes, this might take some time (lines=%1.2e) ...\n\n",$n_of_pfam_parsed_lines);
    if(!$n_of_pfam_parsed_lines)
    {
      die "# EXIT: parsed Pfam output ($pfam_file) seems to be empty, please remove '$newDIR/' and re-run\n";
    }
    construct_Pfam_hash($pfam_file);
    $do_PFAM = $pfam_file;
  }
}

################################################################

# 3) in this section the code deals with homologies based on BLAST similarities      
      
# 3.0) make genome composition report if required
if($do_genome_composition) 
{
  my ($s,$t,$t2,@pangenome,@coregenome,@softcore,$n_of_permutations,$soft_taxa); #$s = sample, $t=taxon to be added, $t2=taxon to compare
  my ($mean,$sd,$data_file,$sort,%previous_sorts,%inparalogues,%homol_registry,@sample,@clusters);
  my @tmptaxa = @taxa;
  my $n_of_taxa = scalar(@tmptaxa);

  if($include_file)
  {
    $NOFSAMPLESREPORT = 1;
    print "\n# genome composition report (samples=1, using sequence order implicit in -I file: $include_file)\n";
  }
  else
  {
    $n_of_permutations = sprintf("%g",factorial($n_of_taxa));
    if($n_of_permutations < $NOFSAMPLESREPORT){ $NOFSAMPLESREPORT = $n_of_permutations; }
    print "\n# genome composition report (samples=$NOFSAMPLESREPORT,permutations=$n_of_permutations,seed=$random_number_generator_seed)\n";
  }

  for($s=0;$s<$NOFSAMPLESREPORT;$s++) # random-sort the list of taxa $NOFSAMPLESREPORT times
  {
    if(!$include_file && $s) # reshuffle until a new permutation is obtained, conserve input order in first sample
    {
      $sort = fisher_yates_shuffle( \@tmptaxa );
      while($previous_sorts{$sort}){ $sort = fisher_yates_shuffle( \@tmptaxa ); }
      $previous_sorts{$sort} = 1;
    }
    push(@{$sample[$s]},@tmptaxa);
  }

  if($do_soft)
  {
    print "# genomic composition parameters: MIN_PERSEQID_HOM=$MIN_PERSEQID_HOM MIN_COVERAGE_HOM=$MIN_COVERAGE_HOM SOFTCOREFRACTION=$SOFTCOREFRACTION ";
  }
  else{ print "# genomic composition parameters: MIN_PERSEQID_HOM=$MIN_PERSEQID_HOM MIN_COVERAGE_HOM=$MIN_COVERAGE_HOM "; }
  print "(set in lib/marfil_homology.pm)\n# genome order:\n";

  for($t=0;$t<$n_of_taxa;$t++){ print "# $t $taxa[$t]\n"; } print "\n";

  if(!-s $marfil_homology::parameter_INP_log || check_different_params('INP',('EVALUE'=>$evalue_cutoff,
        'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'MINBLAST'=>$do_minimal_BDBHs)))
  {
    save_params('INP',('EVALUE'=>$evalue_cutoff,'PI'=>$pi_cutoff,
        'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'MINBLAST'=>$do_minimal_BDBHs));
    $diff_INP_params = 1;
  }
  if(!-s $marfil_homology::parameter_HOM_log || check_different_params('HOM',('EVALUE'=>$evalue_cutoff,
        'PI'=>$MIN_PERSEQID_HOM,'COVER'=>$MIN_COVERAGE_HOM,'MINBLAST'=>$do_minimal_BDBHs)))
  {
    save_params('HOM',('EVALUE'=>$evalue_cutoff,'PI'=>$MIN_PERSEQID_HOM,
        'COVER'=>$MIN_COVERAGE_HOM,'MINBLAST'=>$do_minimal_BDBHs));
    $diff_HOM_params = 1;
  }

  if($runmode eq 'cluster') # precalculate inparalogues/orthologues in cluster
  {
    my ($command,%cluster_PIDs,@to_be_deleted,$clusteroutfile,$clusterlogfile);
    if($doMCL || (!$doPARANOID && !$doMCL && !$doCOG))
    {
      $redo_inp = $reparse_all || $diff_INP_params;
      for(my $j=0;$j<$n_of_taxa;$j++)
      {
        $clusteroutfile = get_makeInparalog_outfilename($taxa[$j]);
        next if(!$redo_inp && -e $clusteroutfile);
        $command = "$ENV{'EXE_INPARA'} -d $newDIR -b $bpo_file -t $taxa[$j] -E $evalue_cutoff ".
          "-S $pi_cutoff -C $pmatch_cutoff -N $neighbor_corr_cutoff -f $redo_inp "; #die "$command\n";
        $clusterlogfile = $clusteroutfile.'.queue';
        push(@to_be_deleted,$clusterlogfile);
        submit_cluster_job($taxa[$j],$command,$clusterlogfile,$newDIR,
          $SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
      }
      if(@to_be_deleted)
      {
        print "# submitting inparalogues jobs to cluster ...\n";
        check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
        print "\n";
        unlink(@to_be_deleted);
        %cluster_PIDs = @to_be_deleted = ();
      }
    }
    $redo_inp = $reparse_all || $diff_HOM_params;
    for($s=0;$s<$NOFSAMPLESREPORT;$s++)
    {
      @tmptaxa = @{$sample[$s]};
      for(my $j=1;$j<$n_of_taxa;$j++)
      {
        for(my $k=$j-1;$k>=0;$k--)
        {
          $clusteroutfile = get_makeHomolog_outfilename($tmptaxa[$j],$tmptaxa[$k]);
          next if(!$redo_inp && -e $clusteroutfile);
          $clusterlogfile = $clusteroutfile.'.queue';
          next if(grep/^$clusterlogfile$/,@to_be_deleted);
          $command = "$ENV{'EXE_HOMOL'} -d $newDIR -b $bpo_file -i $tmptaxa[$j]".
            " -j $tmptaxa[$k] -E $evalue_cutoff -S $MIN_PERSEQID_HOM ".
            " -C $MIN_COVERAGE_HOM -f $redo_inp ";
          push(@to_be_deleted,$clusterlogfile);
          submit_cluster_job('h'.$tmptaxa[$j].'-'.$tmptaxa[$k],$command,$clusterlogfile,
            $newDIR,$SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
        }
      }
    }
    if(@to_be_deleted)
    {
      print "# submitting homologues jobs to cluster ...\n";
      check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
      print "\n";
      unlink(@to_be_deleted);
    }
  }

  if($doMCL)
  {
    if(!-s $parameter_OMCL_log || check_different_params('OMCL',('EVALUE'=>$evalue_cutoff,'FILES'=>$current_files,
          'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'INFLATION'=>$MCLinflation)))
    {
      save_params('OMCL',('EVALUE'=>$evalue_cutoff,'PI'=>$pi_cutoff,'INFLATION'=>$MCLinflation,
          'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'FILES'=>$current_files));
      $diff_OMCL_params = 1;

      #if(!-s $parameter_OMCL_log || check_different_params('OMCL',('INFLATION'=>$MCLinflation))){ $diff_OMCL_params = 1; }
      if(!-s $parameter_OMCL_log || check_different_params('OMCL',('EVALUE'=>$evalue_cutoff,
            'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff)))
      {
        $diff_BDBH_params = 1;
      }
    }
    if($runmode eq 'cluster')
    {
      my ($command,%cluster_PIDs,@to_be_deleted,$clusteroutfile,$clusterlogfile);
      $redo_inp = $reparse_all || $diff_INP_params || $diff_BDBH_params;
      for(my $i=0;$i<$n_of_taxa-1;$i++)
      {
        for(my $j=$i+1;$j<$n_of_taxa;$j++)
        {
          $clusteroutfile = get_makeOrtholog_outfilename('010',$taxa[$i],$taxa[$j]);
          next if(!$redo_inp && -e $clusteroutfile);
          $command = "$ENV{'EXE_ORTHO'} -d $newDIR -b $bpo_file -i $taxa[$i]".
            " -j $taxa[$j] -E $evalue_cutoff -D 0 -S $pi_cutoff -C $pmatch_cutoff".
            " -N $neighbor_corr_cutoff -f $redo_inp -n 1 -l 0 -B 0"; #die "$command\n";
          $clusterlogfile = $clusteroutfile.'.queue';
          push(@to_be_deleted,$clusterlogfile);
          submit_cluster_job($taxa[$i].'-'.$taxa[$j],$command,$clusterlogfile,
            $newDIR,$SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
        }
      }
      if(@to_be_deleted)
      {
        print "# submitting orthologues jobs to cluster ...\n";
        check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
        print "\n";
        unlink(@to_be_deleted);
      }
      if($diff_INP_params){ $diff_INP_params = -1 }
      if($diff_BDBH_params){ $diff_BDBH_params = -1 }
      $diff_OMCL_params = $reparse_all || $diff_OMCL_params;
    }
    else
    {
      $diff_INP_params = $reparse_all || $diff_INP_params;
      $diff_BDBH_params = $reparse_all || $diff_BDBH_params;
      $diff_OMCL_params = $reparse_all || $diff_OMCL_params;
    }

    my ($ref_hash_orths,$ref_hash_orth_taxa) =
      find_OMCL_clusters( $saveRAM,$evalue_cutoff,$pi_cutoff,$pmatch_cutoff,$neighbor_corr_cutoff,
      $MCLinflation,1,$diff_INP_params,$diff_BDBH_params,$diff_OMCL_params,\@taxa );

    %orth_taxa = %{$ref_hash_orth_taxa};
    if(!$do_PFAM){ %orthologues = %{$ref_hash_orths} }
    else
    {
      my $ref_pfam_orths = split_Pfam_clusters( $ref_hash_orths , \%orth_taxa  );
      %orthologues = %{$ref_pfam_orths};
    }
    
    @clusters = keys(%orthologues);

    $orthoMCLdone = 1;
  }
  elsif($doPARANOID) # gets orths + inparalogs mixed
  {
    if(!-e $paranoidDIR){ mkdir($paranoidDIR) }
    my ($ref_hash_orths,$ref_hash_orth_taxa) =
      find_PARANOID_clusters($saveRAM,$bitscore_cutoff,$pmatch_cutoff,$segment_cover_cutoff,
      $paranoidDIR,$reparse_all,@taxa);

    %orth_taxa = %{$ref_hash_orth_taxa};
    if(!$do_PFAM){ %orthologues = %{$ref_hash_orths} }
    else
    {
      my $ref_pfam_orths = split_Pfam_clusters( $ref_hash_orths , \%orth_taxa  );
      %orthologues = %{$ref_pfam_orths};
    }
    
    @clusters = keys(%orthologues);

    $PARANOIDdone = 1;
  }
  elsif($doCOG) # gets orths + inparalogs mixed
  {
    my ($ref_hash_orths,$ref_hash_orth_taxa) =
      find_COGs($segment_cover_cutoff,$evalue_cutoff,$COGmulticluster,$saveRAM,$pwd,$reparse_all,@taxa);

    %orth_taxa = %{$ref_hash_orth_taxa};
    if(!$do_PFAM){ %orthologues = %{$ref_hash_orths} }
    else
    {
      my $ref_pfam_orths = split_Pfam_clusters( $ref_hash_orths , \%orth_taxa  );
      %orthologues = %{$ref_pfam_orths};
    }

    @clusters = keys(%orthologues);

    $COGdone = 1;
  }
  else # BDBH
  {
    if(!-s $parameter_BDBH_log || check_different_params('BDBH',('EVALUE'=>$evalue_cutoff,
          'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'PFAM'=>$do_PFAM)))
    {
      save_params('BDBH',('EVALUE'=>$evalue_cutoff,'PI'=>$pi_cutoff,'PFAM'=>$do_PFAM,
          'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff));

      $diff_BDBH_params = 1;
    }

    if($runmode eq 'cluster')
    {
      my ($command,%cluster_PIDs,@to_be_deleted,$clusteroutfile,$clusterlogfile,$fmask);
      $redo_inp = $reparse_all || $diff_INP_params || $diff_BDBH_params;
      for($s=0;$s<$NOFSAMPLESREPORT;$s++)
      {
        @tmptaxa = @{$sample[$s]};
        for(my $j=1;$j<$n_of_taxa;$j++)
        {
          $fmask = '100'; if($do_PFAM){ $fmask = '101' }
          $clusteroutfile = get_makeOrtholog_outfilename($fmask,$tmptaxa[0],$tmptaxa[$j]);
          next if(!$redo_inp && -e $clusteroutfile);
          $clusterlogfile = $clusteroutfile.'.queue';
          next if(grep/^$clusterlogfile$/,@to_be_deleted);
          $command = "$ENV{'EXE_ORTHO'} -d $newDIR -b $bpo_file -i $tmptaxa[0]".
            " -j $tmptaxa[$j] -E $evalue_cutoff -D $do_PFAM -S $pi_cutoff -C $pmatch_cutoff".
            " -N $neighbor_corr_cutoff -f $redo_inp -n 0 -l 1 -B 1";
          push(@to_be_deleted,$clusterlogfile);
          submit_cluster_job($tmptaxa[0].'-'.$tmptaxa[$j],$command,$clusterlogfile,
            $newDIR,$SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
        }
      }
      if(@to_be_deleted)
      {
        print "# submitting orthologues jobs to cluster ...\n";
        check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
        print "\n";
        unlink(@to_be_deleted);
      }
    }

    $BDBHdone = 1;
  }

  # 3.0.1) sample taxa in random order
  for($s=0;$s<$NOFSAMPLESREPORT;$s++)
  {
    my (%n_of_taxa_in_cluster,%n_of_homs_in_genomes,$sample);
    @tmptaxa = @{$sample[$s]};

    $sample = "## sample $s ($tmptaxa[0] | ";
    for($t=0;$t<$n_of_taxa;$t++)
    {
      $t2=0;
      while($tmptaxa[$t] ne $taxa[$t2]){ $t2++ }
      $sample .= "$t2,";

      # arbitrary trimming
      if(length($sample)>70){ $sample .= '...'; last }
    }
    $sample .= ')';
    print "$sample\n";

    # 3.0.2) calculate pan/core-genome size adding genomes one-by-one
    if($doMCL || $doPARANOID || $doCOG)
    {
      my $n_of_inparalogues = 0;
      foreach $cluster (@clusters)
      {
        foreach $taxon (keys(%{$orth_taxa{$cluster}}))
        {
          next if($taxon ne $tmptaxa[0]);
          $n_of_inparalogues += ($orth_taxa{$cluster}{$taxon}-1);
          last;
        }
      }
      $coregenome[$s][0] = $gindex{$tmptaxa[0]}[2] - $n_of_inparalogues;
      if($do_soft){ $softcore[$s][0] = $coregenome[$s][0] }
      $pangenome[$s][0]  = $coregenome[$s][0];
      print "# adding $tmptaxa[0]: core=$coregenome[$s][0] pan=$pangenome[$s][0]\n";
    }
    else # BDBH
    {
      if(!$LSE_registry{$tmptaxa[0]} &&
        ($runmode eq 'cluster' && !-e get_makeInparalog_outfilename($tmptaxa[0])))
      {
        $redo_inp = $reparse_all || $diff_INP_params;
        $LSE_registry{$tmptaxa[0]} = 1;
      }
      elsif($runmode eq 'local')
      {
        $redo_inp = $reparse_all || $diff_INP_params;
      }
      else{ $redo_inp = 0 }

      print "# clustering inparalogues in $tmptaxa[0]\n";
      my($ref_inparalogues) = makeInparalog($saveRAM,$tmptaxa[0],$evalue_cutoff,
        $pi_cutoff,$pmatch_cutoff,$neighbor_corr_cutoff,1,$redo_inp);
      $LSE_reference = cluster_lineage_expansions($ref_inparalogues);
      
      # calculate initial core & pan size
      my $initial_core_size = 0; 
      foreach $gene ($gindex{$tmptaxa[0]}[0] .. $gindex{$tmptaxa[0]}[1])
      {
        next if( $LSE_reference->{$gene} ); # inparalogues
 
        $initial_core_size++;
      }
           
      $coregenome[$s][0] = $initial_core_size;
      
      if($coregenome[$s][0] < 0){ $coregenome[$s][0] = 0 }
      
      $pangenome[$s][0]  = $coregenome[$s][0];
      print "# adding $tmptaxa[0]: core=$coregenome[$s][0] pan=$pangenome[$s][0]\n";
    }

    for($t=1;$t<$n_of_taxa;$t++)
    {
      $coregenome[$s][$t] = 0;
      $pangenome[$s][$t] = $pangenome[$s][$t-1];

      # core genome
      if($doMCL || $doPARANOID || $doCOG)
      {
        CLUSTER: foreach $cluster (@clusters)
        {
          # potential core clusters must contain sequences from reference taxon $tmptaxa[0]
          # this check is done only once ($t=1)
          next if($t == 1 && !$do_soft && !$orth_taxa{$cluster}{$tmptaxa[0]});
          
          foreach $taxon (keys(%{$orth_taxa{$cluster}}))
          {
            if($taxon eq $tmptaxa[$t])
            {
              $n_of_taxa_in_cluster{$cluster}++; # taxa added starting from $t=1
              
              if($orth_taxa{$cluster}{$tmptaxa[0]} && $n_of_taxa_in_cluster{$cluster} == $t)
              {
                $coregenome[$s][$t]++;  # update core totals
                if($do_soft){ $softcore[$s][$t]++ }
              }
              elsif($do_soft) 
              {
                $soft_taxa = $n_of_taxa_in_cluster{$cluster} || 0;
                if($orth_taxa{$cluster}{$tmptaxa[0]}){ $soft_taxa++ }
                
                if($soft_taxa >= int(($t+1)*$SOFTCOREFRACTION))
                {
                  $softcore[$s][$t]++;  
                }
              }
                     
              next CLUSTER;
            }
          }
        }
        #print "# adding $tmptaxa[$t]: core=$coregenome[$s][$t] pan=$pangenome[$s][$t]\n";
      }
      else # BDBH core: orthologues/pairs among two genomes (transitivity: if(orth(0,1) && orth(0,2)) then orth(1,2)
      {
        my $ref_orths;
        $label = $tmptaxa[0].' '.$tmptaxa[$t];

        if(!$LSE_registry{$tmptaxa[$t]} &&
          ($runmode eq 'cluster' && !-e get_makeInparalog_outfilename($tmptaxa[$t])))
        {
          $redo_inp = $reparse_all || $diff_INP_params;
          $LSE_registry{$tmptaxa[$t]} = 1;
        }
        elsif($runmode eq 'local')
        {
          $redo_inp = $reparse_all || $diff_INP_params;
        }
        else{ $redo_inp = 0 }

        my $fmask = '100'; if($do_PFAM){ $fmask = '101' }
        if(!$orth_registry{$label} && ($runmode eq 'cluster' &&
            !-e get_makeOrtholog_outfilename($fmask,$tmptaxa[0],$tmptaxa[$t])))
        {
          $redo_orth = $diff_BDBH_params;
          $orth_registry{$label} = 1;
        }
        elsif($runmode eq 'local')
        {
          $redo_orth = $diff_BDBH_params;
        }
        else{ $redo_orth = 0 }

        print "# clustering inparalogues in $tmptaxa[$t]\n";
        my($ref_inparalogues) = makeInparalog($saveRAM,$tmptaxa[$t],$evalue_cutoff,$pi_cutoff,
          $pmatch_cutoff,$neighbor_corr_cutoff,1,$redo_inp);
        $LSE_t = cluster_lineage_expansions($ref_inparalogues);
        %inparalogues = %$LSE_t;

        print "# finding BDBHs between $tmptaxa[0] and $tmptaxa[$t]\n";
        @{$ref_orths} = makeOrtholog($saveRAM,$tmptaxa[0],$tmptaxa[$t],
          1,$evalue_cutoff,$pi_cutoff,
          $pmatch_cutoff,$neighbor_corr_cutoff,0,$do_PFAM,$redo_orth,$LSE_reference,$LSE_t);

        foreach $gene ($gindex{$tmptaxa[0]}[0] .. $gindex{$tmptaxa[0]}[1])
        {
          if($ref_orths->[0]->{$gene})
          {
            $n_of_homs_in_genomes{$gene}++;
            if($n_of_homs_in_genomes{$gene} == $t)
            {
              $coregenome[$s][$t]++; # update core
            }
          }
        }
      }

      # pan genome (unique genes) : those without hits in last added genome when compared to all previous
      for($t2=$t-1;$t2>=0;$t2--)
      {
        $label = $tmptaxa[$t].' '.$tmptaxa[$t2];
        if(!$homol_registry{$label} && ($runmode eq 'cluster' &&
            !-e get_makeHomolog_outfilename($tmptaxa[$t],$tmptaxa[$t2])))
        {
          $redo_orth = $diff_HOM_params;
          $homol_registry{$label} = 1;
        }
        elsif($runmode eq 'local')
        {
          $redo_orth = $diff_HOM_params;
        }
        else{ $redo_orth = 0 }

        print "# finding homologs between $tmptaxa[$t] and $tmptaxa[$t2]\n";
        my $ref_homol = makeHomolog($saveRAM,$tmptaxa[$t],$tmptaxa[$t2],$evalue_cutoff,
          $MIN_PERSEQID_HOM,$MIN_COVERAGE_HOM,$redo_orth);

        foreach $gene ($gindex{$tmptaxa[$t]}[0] .. $gindex{$tmptaxa[$t]}[1])
        {
          if($ref_homol->{$gene}){ $n_of_homs_in_genomes{$gene}++; }
        }
      }

      # label inparalogues in OMCL,PRND,COGS clusters to avoid over-estimating pangenome
      if($doMCL || $doPARANOID || $doCOG)
      {
        %inparalogues = ();
        foreach $cluster (@clusters)
        {
          # skip clusters with <2 $t sequences
          next if(!$orth_taxa{$cluster}{$tmptaxa[$t]} ||
            $orth_taxa{$cluster}{$tmptaxa[$t]} < 2);
            
          foreach $gene (@{$orthologues{$cluster}})
          {
            next if($gindex2[$gene] ne $tmptaxa[$t]);
            $inparalogues{$gene} = 1;
          }
        }
      }

      # update pan total
      foreach $gene ($gindex{$tmptaxa[$t]}[0] .. $gindex{$tmptaxa[$t]}[1])
      {
        next if($n_of_homs_in_genomes{$gene} || $inparalogues{$gene});
        
        $pangenome[$s][$t]++;
      }

      print "# adding $tmptaxa[$t]: core=$coregenome[$s][$t] pan=$pangenome[$s][$t]\n"; 
    }
  }

  # 3.0.3) print pan-genome composition stats
  $data_file = $newDIR ."/pan_genome".$pancore_mask.".tab";
  print "\n# pan-genome (number of genes, can be plotted with plot_pancore_matrix.pl)\n# file=".
    short_path($data_file,$pwd)."\n";
  print "genomes\tmean\tstddev\t|\tsamples\n";
  for($t=0;$t<$n_of_taxa;$t++)
  {
    my @data;
    for($s=0;$s<$NOFSAMPLESREPORT;$s++){ push(@data,$pangenome[$s][$t]) }
    $mean = sprintf("%1.0f",calc_mean(\@data));
    $sd = sprintf("%1.0f",calc_std_deviation(\@data));
    print "$t\t$mean\t$sd\t|\t";
    for($s=0;$s<$NOFSAMPLESREPORT;$s++){ print "$pangenome[$s][$t]\t"; } print "\n";
  }

  # 3.0.4) create input file for pan-genome composition boxplot
  open(BOXDATA,">$data_file") || die "# EXIT: cannot create $data_file\n";
  for($t=0;$t<$n_of_taxa;$t++)
  {
    $label = 'g'.($t+1);
    print BOXDATA "$label\t";
  } print BOXDATA "\n";

  for($s=0;$s<$NOFSAMPLESREPORT;$s++)
  {
    for($t=0;$t<$n_of_taxa;$t++){ print BOXDATA "$pangenome[$s][$t]\t";}
    print BOXDATA "\n";
  }
  close(BOXDATA);

  # 3.0.5) print core-genome composition stats
  $data_file = $newDIR ."/core_genome".$pancore_mask.".tab";
  print "\n# core-genome (number of genes, can be plotted with plot_pancore_matrix.pl)\n# file=".
    short_path($data_file,$pwd)."\n";
  print "genomes\tmean\tstddev\t|\tsamples\n";
  for($t=0;$t<$n_of_taxa;$t++)
  {
    my @data;
    for($s=0;$s<$NOFSAMPLESREPORT;$s++){ push(@data,$coregenome[$s][$t]) }
    $mean = sprintf("%1.0f",calc_mean(\@data));
    $sd = sprintf("%1.0f",calc_std_deviation(\@data));
    print "$t\t$mean\t$sd\t|\t";
    for($s=0;$s<$NOFSAMPLESREPORT;$s++){ print "$coregenome[$s][$t]\t"; } print "\n";
  }

  # 3.0.6) create input file for core-genome composition boxplot
  open(BOXDATA,">$data_file") || die "# EXIT : cannot create $data_file\n";
  for($t=0;$t<$n_of_taxa;$t++)
  {
    $label = 'g'.($t+1);
    print BOXDATA "$label\t";
  } print BOXDATA "\n";

  for($s=0;$s<$NOFSAMPLESREPORT;$s++)
  {
    for($t=0;$t<$n_of_taxa;$t++){ print BOXDATA "$coregenome[$s][$t]\t";}
    print BOXDATA "\n";
  }
  close(BOXDATA);
  
  if($do_soft)
  {
    # 3.0.7) print soft core-genome composition stats
    $data_file = $newDIR ."/soft-core_genome".$pancore_mask.".tab";
    print "\n# soft core-genome (number of genes, can be plotted with plot_pancore_matrix.pl)\n# file=".
      short_path($data_file,$pwd)."\n";
    print "genomes\tmean\tstddev\t|\tsamples\n";
    for($t=0;$t<$n_of_taxa;$t++)
    {
      my @data;
      for($s=0;$s<$NOFSAMPLESREPORT;$s++){ push(@data,$softcore[$s][$t]) }
      $mean = sprintf("%1.0f",calc_mean(\@data));
      $sd = sprintf("%1.0f",calc_std_deviation(\@data));
      print "$t\t$mean\t$sd\t|\t";
      for($s=0;$s<$NOFSAMPLESREPORT;$s++){ print "$softcore[$s][$t]\t"; } print "\n";
    }
  
    # 3.0.8) create input file for soft core-genome composition boxplot
    open(BOXDATA,">$data_file") || die "# EXIT : cannot create $data_file\n";
    for($t=0;$t<$n_of_taxa;$t++)
    {
      $label = 'g'.($t+1);
      print BOXDATA "$label\t";
    } print BOXDATA "\n";
  
    for($s=0;$s<$NOFSAMPLESREPORT;$s++)
    {
      for($t=0;$t<$n_of_taxa;$t++){ print BOXDATA "$softcore[$s][$t]\t";}
      print BOXDATA "\n";
    }
    close(BOXDATA);
  }

  # add missing singleton clusters for later use, they are not needed for -c
  if($min_cluster_size == 0 && ($doMCL || $doPARANOID || $doCOG))
  {
    add_unmatched_singletons(\%orthologues,\%orth_taxa,\@taxa);
  }
}

# 3.1) identify orthologues and inparalogues using reference/OMCL/PRND paths #############################
print "\n# clustering orthologous sequences\n";

if(!-s $marfil_homology::parameter_INP_log || check_different_params('INP',('EVALUE'=>$evalue_cutoff,
      'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'MINBLAST'=>$do_minimal_BDBHs)))
{
  save_params('INP',('EVALUE'=>$evalue_cutoff,'PI'=>$pi_cutoff,
      'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'MINBLAST'=>$do_minimal_BDBHs));
  $diff_INP_params = 1;
}

#for(my $i=0;$i<$n_of_taxa;$i++){
#  next if(!$cluster_ids{$taxa[$i]});
#  print "\n# clustering inparalogues in $taxa[$i] (cluster members)\n";
#  makeClusterInparalog($saveRAM,$taxa[$i],\%cluster_ids,$reparse_all);          
#}

if($runmode eq 'cluster')
{
  my ($command,%cluster_PIDs,@to_be_deleted,$clusteroutfile,$clusterlogfile);
  if(($doMCL && !$orthoMCLdone) || (!$doPARANOID && !$doMCL && !$doCOG))
  {
    $redo_inp = ($reparse_all || $diff_INP_params) && !$BDBHdone;
    for(my $j=0;$j<$n_of_taxa;$j++)
    {
      $clusteroutfile = get_makeInparalog_outfilename($taxa[$j]);
      next if(!$redo_inp && -e $clusteroutfile);
      $command = "$ENV{'EXE_INPARA'} -d $newDIR -b $bpo_file -t $taxa[$j] -E $evalue_cutoff ".
        "-S $pi_cutoff -C $pmatch_cutoff -N $neighbor_corr_cutoff -f $redo_inp ";
      $clusterlogfile = $clusteroutfile.'.queue';
      push(@to_be_deleted,$clusterlogfile);
      submit_cluster_job($taxa[$j],$command,$clusterlogfile,$newDIR,
        $SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
    }
    if(@to_be_deleted)
    {
      print "\n# submitting inparalogues jobs to cluster ...\n";
      check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
      unlink(@to_be_deleted);
      %cluster_PIDs = @to_be_deleted = ();
    }
  }
}

if($doMCL)
{
  if(!$orthoMCLdone)
  {
    # find out which previous results, if any, can be recyled
    if(!-s $parameter_OMCL_log || check_different_params('OMCL',('EVALUE'=>$evalue_cutoff,'FILES'=>$current_files,
          'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'INFLATION'=>$MCLinflation)))
    {
      save_params('OMCL',('EVALUE'=>$evalue_cutoff,'PI'=>$pi_cutoff,'INFLATION'=>$MCLinflation,
          'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'FILES'=>$current_files));
      $diff_OMCL_params = 1;

      #if(!-s $parameter_OMCL_log || check_different_params('OMCL',('INFLATION'=>$MCLinflation))){ $diff_OMCL_params = 1; }
      if(!-s $parameter_OMCL_log || check_different_params('OMCL',('EVALUE'=>$evalue_cutoff,
            'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff)))
      {
        $diff_BDBH_params = 1;
      }
    }

    if($runmode eq 'cluster')
    {
      my ($command,%cluster_PIDs,@to_be_deleted,$clusteroutfile,$clusterlogfile);
      $redo_inp = $reparse_all || $diff_INP_params || $diff_BDBH_params;
      for(my $i=0;$i<$n_of_taxa-1;$i++)
      {
        for(my $j=$i+1;$j<$n_of_taxa;$j++)
        {
          $clusteroutfile = get_makeOrtholog_outfilename('010',$taxa[$i],$taxa[$j]);
          next if(!$redo_inp && -e $clusteroutfile);
          $command = "$ENV{'EXE_ORTHO'} -d $newDIR -b $bpo_file -i $taxa[$i]".
            " -j $taxa[$j] -E $evalue_cutoff -D 0 -S $pi_cutoff -C $pmatch_cutoff".
            " -N $neighbor_corr_cutoff -f $redo_inp -n 1 -l 0 -B 0"; #die "$command\n";
          $clusterlogfile = $clusteroutfile.'.queue';
          push(@to_be_deleted,$clusterlogfile);
          submit_cluster_job($taxa[$i].'-'.$taxa[$j],$command,$clusterlogfile,
            $newDIR,$SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
        }
      }
      if(@to_be_deleted)
      {
        print "\n# submitting orthologues jobs to cluster ...\n";
        check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
        unlink(@to_be_deleted);
      }
      if($diff_INP_params){ $diff_INP_params = -1 }
      if($diff_BDBH_params){ $diff_BDBH_params = -1 }
      $diff_OMCL_params = $reparse_all || $diff_OMCL_params;
    }
    else
    {
      $diff_INP_params = $reparse_all || $diff_INP_params;
      $diff_BDBH_params = $reparse_all || $diff_BDBH_params;
      $diff_OMCL_params = $reparse_all || $diff_OMCL_params;
    }

    my ($ref_hash_orths,$ref_hash_orth_taxa) =
      find_OMCL_clusters( $saveRAM,$evalue_cutoff,$pi_cutoff,$pmatch_cutoff,$neighbor_corr_cutoff,
      $MCLinflation,1,$diff_INP_params,$diff_BDBH_params,$diff_OMCL_params,\@taxa );

    if($min_cluster_size == 0)
    {
      add_unmatched_singletons($ref_hash_orths,$ref_hash_orth_taxa,\@taxa);
    }

    %orth_taxa = %{$ref_hash_orth_taxa};
    if(!$do_PFAM){ %orthologues = %{$ref_hash_orths} }
    else
    {
      my $ref_pfam_orths = split_Pfam_clusters( $ref_hash_orths , \%orth_taxa );
      %orthologues = %{$ref_pfam_orths};
    }
  }
}
elsif($doPARANOID)
{
  if(!$PARANOIDdone)
  {
    if(!-e $paranoidDIR){ mkdir($paranoidDIR) }
    my ($ref_hash_orths,$ref_hash_orth_taxa) =
      find_PARANOID_clusters($saveRAM,$bitscore_cutoff,$pmatch_cutoff,
      $segment_cover_cutoff,$paranoidDIR,$reparse_all,@taxa);

    if($min_cluster_size == 0)
    {
      add_unmatched_singletons($ref_hash_orths,$ref_hash_orth_taxa,\@taxa);
    }

    %orth_taxa = %{$ref_hash_orth_taxa};
    if(!$do_PFAM){ %orthologues = %{$ref_hash_orths} }
    else
    {
      my $ref_pfam_orths = split_Pfam_clusters( $ref_hash_orths , \%orth_taxa  );
      %orthologues = %{$ref_pfam_orths};
    }
  }
}
elsif($doCOG)
{
  if(!$COGdone)
  {
    my ($ref_hash_orths,$ref_hash_orth_taxa) =
      find_COGs($segment_cover_cutoff,$evalue_cutoff,$COGmulticluster,$saveRAM,$pwd,$reparse_all,@taxa);

    if($min_cluster_size == 0)
    {
      add_unmatched_singletons($ref_hash_orths,$ref_hash_orth_taxa,\@taxa);
    }

    %orth_taxa = %{$ref_hash_orth_taxa};
    if(!$do_PFAM){ %orthologues = %{$ref_hash_orths} }
    else
    {
      my $ref_pfam_orths = split_Pfam_clusters( $ref_hash_orths , \%orth_taxa , 1 );
      %orthologues = %{$ref_pfam_orths};
    }
  }
}
else # BDBH
{
  if(!-s $parameter_BDBH_log || check_different_params('BDBH',('EVALUE'=>$evalue_cutoff,
        'PI'=>$pi_cutoff,'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff,'PFAM'=>$do_PFAM)))
  {
    save_params('BDBH',('EVALUE'=>$evalue_cutoff,'PI'=>$pi_cutoff,'PFAM'=>$do_PFAM,
        'COVER'=>$pmatch_cutoff,'NCORR'=>$neighbor_corr_cutoff));
    $diff_BDBH_params = 1;
  }#print "$diff_INP_params $diff_BDBH_params\n";

  if($runmode eq 'cluster') # precalculate orthologues in cluster
  {
    my ($command,%cluster_PIDs,@to_be_deleted,$clusteroutfile,$clusterlogfile,$fmask);
    $redo_inp = ($reparse_all || $diff_INP_params || $diff_BDBH_params) && !$BDBHdone;
    for(my $j=0;$j<$n_of_taxa;$j++)
    {
      next if($j == $reference_proteome);
      $fmask = '100'; if($do_PFAM){ $fmask = '101' }
      $clusteroutfile = get_makeOrtholog_outfilename($fmask,$taxa[$reference_proteome],$taxa[$j]);
      next if(!$redo_inp && -e $clusteroutfile);
      $command = "$ENV{'EXE_ORTHO'} -d $newDIR -b $bpo_file -i $taxa[$reference_proteome]".
        " -j $taxa[$j] -E $evalue_cutoff -D $do_PFAM -S $pi_cutoff -C $pmatch_cutoff".
        " -N $neighbor_corr_cutoff -f $redo_inp -n 0 -l 1 -B 1"; #die "$command\n";
      $clusterlogfile = $clusteroutfile.'.queue';
      push(@to_be_deleted,$clusterlogfile);
      submit_cluster_job($taxa[$reference_proteome].'-'.$taxa[$j],$command,$clusterlogfile,
        $newDIR,$SGEPATH,$QUEUESETTINGS,$QUEUEWAIT,\%cluster_PIDs);
    }
    if(@to_be_deleted)
    {
      print "\n# submitting orthologues jobs to cluster ...\n";
      check_cluster_jobs($newDIR,$SGEPATH,$QUEUEWAIT,$WAITTIME,\%cluster_PIDs);
      unlink(@to_be_deleted);
    }
  }

  # 3.1.1) cluster lineage-specificas expansions (LSEs, inparalogs) in reference proteome
  print "\n# clustering inparalogues in $taxa[$reference_proteome] (reference)\n";
  if(!$LSE_registry{$taxa[$reference_proteome]} &&
    ($runmode eq 'cluster' && !-e get_makeInparalog_outfilename($taxa[$reference_proteome])))
  {
    $redo_inp = $reparse_all || $diff_INP_params;
    $LSE_registry{$taxa[$reference_proteome]} = 1;
  }
  elsif($runmode eq 'local')
  {
    $redo_inp = $reparse_all || $diff_INP_params;
  }
  else{ $redo_inp = 0 }

  my($rhash_inparalogues_ref) = makeInparalog($saveRAM,$taxa[$reference_proteome],$evalue_cutoff,
    $pi_cutoff,$pmatch_cutoff,$neighbor_corr_cutoff,1,$redo_inp);
  $LSE_reference = cluster_lineage_expansions($rhash_inparalogues_ref);

  for(my $j=0;$j<$n_of_taxa;$j++)
  {
    next if($j == $reference_proteome);
    $label = $taxa[$reference_proteome].' '.$taxa[$j];

    # 3.1.2) find LSEs within $taxa[$j]
    print "\n# clustering inparalogues in $taxa[$j]\n";
    if(!$LSE_registry{$taxa[$j]} &&
      ($runmode eq 'cluster' && !-e get_makeInparalog_outfilename($taxa[$j])))
    {
      $redo_inp = $reparse_all || $diff_INP_params;
      $LSE_registry{$taxa[$j]} = 1;
    }
    elsif($runmode eq 'local')
    {
      $redo_inp = $reparse_all || $diff_INP_params;
    }
    else{ $redo_inp = 0 }

    my($rhash_inparalogues_j) = makeInparalog($saveRAM,$taxa[$j],$evalue_cutoff,$pi_cutoff,
      $pmatch_cutoff,$neighbor_corr_cutoff,1,$redo_inp);
    $LSE_t = cluster_lineage_expansions($rhash_inparalogues_j);

    print "\n# finding BDBHs between $taxa[$reference_proteome] and $taxa[$j]\n";
    my $fmask = '100'; if($do_PFAM){ $fmask = '101' }
    if(!$orth_registry{$label} &&
      ($runmode eq 'cluster' &&
        !-e get_makeOrtholog_outfilename($fmask,$taxa[$reference_proteome],$taxa[$j])))
    {
      $redo_orth = $reparse_all || $diff_INP_params || $diff_BDBH_params;
      $orth_registry{$label} = 1;
    }
    elsif($runmode eq 'local')
    {
      $redo_orth = $reparse_all || $diff_INP_params || $diff_BDBH_params;
    }
    else{ $redo_orth = 0 }

    my($orth_table_ref_j) = makeOrtholog($saveRAM,$taxa[$reference_proteome],$taxa[$j],1,$evalue_cutoff,$pi_cutoff,
      $pmatch_cutoff,$neighbor_corr_cutoff,0,$do_PFAM,$redo_orth,$LSE_reference,$LSE_t);

    foreach $gene (keys(%$orth_table_ref_j)) # each seed $gene in the reference genome
    {

     #next if($exclude_inparalogues && grep(/^$gene$/,values(%$LSE_reference)));
      foreach $orth (@{$orth_table_ref_j->{$gene}})
      {

        #next if($exclude_inparalogues && grep(/^$orth$/,values(%$LSE_t)));
        push(@{$orthologues{$gene}},$orth); #print "# $gene $orth\n" if($gene == 77); # debug

        # 3.1.3) add inparalogs from j-th proteome
        $orth_taxa{$gene}{$taxa[$j]} = 1;
        foreach $inpara (keys(%$LSE_t))
        {
          next if($LSE_t->{$inpara} ne $orth);
          push(@{$orthologues{$gene}},$inpara); #print "## j $inpara ($orth)\n";
          $orth_taxa{$gene}{$taxa[$j]}++;
        }
      }
    }
  }

  # 3.1.4) add inparalogs from reference proteome to any existing cluster
  foreach $gene (keys(%orthologues))
  {
    $orth_taxa{$gene}{$taxa[$reference_proteome]} = 1; # accounts for reference proteome seed
    foreach $inpara (keys(%$LSE_reference))
    {
      next if($LSE_reference->{$inpara} ne $gene);
      unshift(@{$orthologues{$gene}},$inpara);
      $orth_taxa{$gene}{$taxa[$reference_proteome]}++;
    }
  }
}

# 4) print clusters of sequences to appropriate directories ################

if($do_features)
{
  print "\n# looking for valid sequence clusters (n_of_taxa=$min_cluster_size)...\n";
  print "# considered GenBank features: $do_features\n\n";
}
else
{
  print "\n# looking for valid ORF clusters (n_of_taxa=$min_cluster_size)...\n\n";
}

$FASTAresultsDIR   = $newDIR ."/".$output_mask;
$cluster_list_file = $FASTAresultsDIR .".cluster_list";

if($do_ANIb_matrix){ $ANIb_matrix_file = $FASTAresultsDIR.'Avg_identity.tab'; }

open(CLUSTERLIST,">$cluster_list_file") || die "# EXIT: cannot create $cluster_list_file\n";

$n_of_clusters = 0;
GENE: foreach $gene (sort {$a<=>$b} (keys(%orthologues)))
{
  $n_of_taxa_cluster = scalar(keys(%{$orth_taxa{$gene}}));
  next if($n_of_taxa_cluster < $min_cluster_size);

  if($exclude_inparalogues)
  {
    foreach $taxon (keys(%{$orth_taxa{$gene}}))
    {
      next GENE if($orth_taxa{$gene}{$taxon} > 1);
    }
  }

  my (@taxon_names,$cluster_name,$header,%aligned_coords,%cluster_taxa,$ref_hash_short_orthologues);

  # 4.1) check length of sequences in clusters if requested
  if($filter_by_length)
  {
    $n_of_similar_length_orthologues=$n_of_similar_length_paralogues=0;

    my @orths = ($gene,@{$orthologues{$gene}});

    ($n_of_similar_length_orthologues,$ref_hash_short_orthologues) = same_length(\%seq_length,\@orths,$filter_by_length);

    next if($n_of_similar_length_orthologues < $min_cluster_size);
  }

  # 4.2) create output directory if necessary
  if(!-e $FASTAresultsDIR){ mkdir($FASTAresultsDIR); }

  # 4.3) print cluster file
  $cluster_size = $dna_cluster_size = 0;
  $prot_cluster = $dna_cluster = '';
  $annot = $pfam_annot = '';
  $dna_clusterfile = 'void';

  # 4.3.0) set cluster name according to reference genome if possible
  if($gindex2[$gene] ne $taxa[$reference_proteome])
  {
    foreach $orth (@{$orthologues{$gene}})
    {
      next if($gindex2[$orth] ne $taxa[$reference_proteome]);
      $cluster_name = extract_gene_name($sequence_data[$orth]);
      if(!$cluster_name){ $cluster_name = (split(/\s+/,$sequence_data[$orth]))[0] }
      $generef = $orth;
      last;
    }

    if(!$cluster_name)
    {
      $cluster_name = extract_gene_name($sequence_data[$gene]);
      if(!$cluster_name){ $cluster_name = (split(/\s+/,$sequence_data[$gene]))[0] }
      $generef = $gene;
    }
  }
  else
  {
    $cluster_name = extract_gene_name($sequence_data[$gene]);
    if(!$cluster_name){ $cluster_name = (split(/\s+/,$sequence_data[$gene]))[0] }
    $generef = $gene;
  }

  $cluster_name =~ s/ /_/g;
  $cluster_name =~ s/[\s|\(|\)|\*|;|\|\[|\]|\/|:|,]/-/g;
  $cluster_name = $generef.'_'.$cluster_name;
  if($names{$cluster_name})
  {
    # avoid unlikely but possible name repeats (COGS)
    $cluster_name = sprintf("%s.%d",$cluster_name,++$names{$cluster_name});
  }
  else{ $names{$cluster_name}++ }

  if(!$saveRAM){ %aligned_coords = find_local_alignment_coords( $gene,@{$orthologues{$gene}} ); }
  if($do_PFAM){ $pfam_annot = $pfam_hash{$gene} || ''; }

  # 4.3.1) print reference gene/protein
  push(@taxon_names,$gindex2[$gene]); # global set in phyTools:constructAllFasta
  $header = $sequence_data[$gene];
  if(!$saveRAM)
  {
    chomp($header);
    $header .= " | aligned:$aligned_coords{$gene}{'first'}-$aligned_coords{$gene}{'last'} ".
      "($aligned_coords{$gene}{'length'})";
  }
  $prot_cluster .= ">$header\n".$sequence_prot[$gene]."\n";
  $cluster_size++;

  if($sequence_dna[$gene])
  {
    $dna_cluster .= ">$sequence_data[$gene]\n".$sequence_dna[$gene]."\n";
    $dna_cluster_size++;
  }
  
  if($total_clustersOK)
  { 
    if($header =~ /\[(.*)?\]/){ $cluster_taxa{$1}++ }
    
    # check for redundant pre-clustered sequences
    if($cluster_ids{$gene})
    {
      foreach $clgene (@{$cluster_ids{$gene}})
      {
        $header = $sequence_data[$clgene]; chomp($header);
        $prot_cluster .= ">$header (pre-clustered)\n".$sequence_prot[$clgene]."\n";
        $cluster_size++;
        
        if($sequence_dna[$clgene])
        {
          $dna_cluster .= ">$sequence_data[$clgene] (pre-clustered)\n".$sequence_dna[$clgene]."\n";
          $dna_cluster_size++;
        }
        if($header =~ /\[(.*)?\]/){ $cluster_taxa{$1}++ }
      }
    }
  }
  elsif($do_intergenic)
  {
    my $gi = (split(/\s+/,$header))[0];  #$gi =~ s/>//;
    push(@{$GIclusters{$cluster_name}},$gi);
    push(@{$idclusters{$cluster_name}},$gene);
  }
  
  # 4.3.2) print orthologues
  foreach $orth (@{$orthologues{$gene}})
  {
    next if($ref_hash_short_orthologues->{$orth});
    push(@taxon_names,$gindex2[$orth]); 

    $header = $sequence_data[$orth];
    if(!$saveRAM)
    {
      chomp($header);
      $header .= " | aligned:$aligned_coords{$orth}{'first'}-$aligned_coords{$orth}{'last'} ".
        "($aligned_coords{$orth}{'length'})";
    }
    $prot_cluster .= ">$header\n".$sequence_prot[$orth]."\n";
    $cluster_size++;

    if($sequence_dna[$orth])
    {
      $dna_cluster .= ">$sequence_data[$orth]\n".$sequence_dna[$orth]."\n";
      $dna_cluster_size++;
    }
    
    if($total_clustersOK)
    { 
      if($header =~ /\[(.*)?\]/){ $cluster_taxa{$1}++ }
    
      # check for redundant pre-clustered sequences
      if($cluster_ids{$orth})
      {
        foreach $clorth (@{$cluster_ids{$orth}})
        {
          $header = $sequence_data[$clorth]; chomp($header);
          $prot_cluster .= ">$header (pre-clustered)\n".$sequence_prot[$clorth]."\n";
          $cluster_size++;
          
          if($sequence_dna[$clorth])
          {
            $dna_cluster .= ">$sequence_data[$clorth] (pre-clustered)\n".$sequence_dna[$clorth]."\n";
            $dna_cluster_size++;
          }
          if($header =~ /\[(.*)?\]/){ $cluster_taxa{$1}++ }
        }
      }
    }
    elsif($do_intergenic)
    {
      my $gi = (split(/\s+/,$header))[0]; #$gi =~ s/>//;
      push(@{$GIclusters{$cluster_name}},$gi);
      push(@{$idclusters{$cluster_name}},$orth);
    }
  }

  $clusterfile = $FASTAresultsDIR ."/$cluster_name.faa";
  if($do_features){ $clusterfile = $FASTAresultsDIR ."/$cluster_name.fna"; }
  open(CLUSTER,">$clusterfile") || die "# EXIT : cannot create $clusterfile\n";
  print CLUSTER $prot_cluster;
  close(CLUSTER);

  # 4.3.4) check dna cluster and print it
  # Error with CONTIG feature in GenBank files (Jul2014): !0 && 68 == 74 
  if(!$do_features && $dna_cluster_size == $cluster_size)
  {
    $dna_clusterfile = $FASTAresultsDIR ."/$cluster_name.fna";
    open(DNACLUSTER,">$dna_clusterfile") || die "# EXIT : cannot create $dna_clusterfile\n";
    print DNACLUSTER $dna_cluster;
    close(DNACLUSTER);
  }

  # 4.3.5) print cluster info
  $annot = "cluster $cluster_name size=$cluster_size taxa=$n_of_taxa_cluster ";
  if($total_clustersOK && %cluster_taxa){ $annot .= sprintf("real_taxa=%d ",scalar(keys(%cluster_taxa))); }  
  if($pfam_annot){ $annot .= "Pfam=$pfam_annot "; }
  if($do_features){ $annot .= "file: void dnafile: $cluster_name.fna\n"; }
  else
  {
    $annot .= "file: $cluster_name.faa ";
    if($dna_clusterfile ne 'void'){ $dna_clusterfile = $cluster_name.'.fna'; }
    $annot .= "dnafile: $dna_clusterfile\n";
  }

  print $annot if($PRINTCLUSTERSSCREEN);

  print CLUSTERLIST $annot;
  foreach $taxon (@taxon_names){ print CLUSTERLIST ": $taxon\n"; }

  # 4.3.6) add data to ANIb matrix if required
  if($do_ANIb_matrix)
  {
    my @genes = sort {$a<=>$b} ($gene,@{$orthologues{$gene}});
    foreach $orth (0 .. $#genes)
    {
      foreach my $orth2 ($orth+1 .. $#genes)
      {
        my @blast_data = blastqueryab($genes[$orth],$genes[$orth2]);
        next if(!$blast_data[3]);
        push(@{$ANIb_matrix{$taxon_names[$orth]}{$taxon_names[$orth2]}},$blast_data[3]);
      }
    }
  }

  $n_of_clusters++;
}

if(!$n_of_clusters){ print "\n# number_of_clusters = $n_of_clusters\n"; }
else
{
  print "\n# number_of_clusters = $n_of_clusters\n# cluster_list = ".short_path($cluster_list_file,$pwd)."\n".
    "# cluster_directory = ".short_path($FASTAresultsDIR,$pwd)."\n";
}

close(CLUSTERLIST);




if($do_ANIb_matrix)
{
  open(ANIBMATRIX,">$ANIb_matrix_file") || die "# EXIT: cannot create $ANIb_matrix_file\n";

  print ANIBMATRIX "genomes";
  for($taxon=0;$taxon<scalar(@taxa);$taxon++)
  {
    print ANIBMATRIX "\t$taxa[$taxon]";
  } print ANIBMATRIX "\n";

  for($taxon=0;$taxon<scalar(@taxa);$taxon++)
  {
    print ANIBMATRIX "$taxa[$taxon]";
    for(my $taxon2=0;$taxon2<scalar(@taxa);$taxon2++)
    {
      if($taxon == $taxon2){ print ANIBMATRIX "\t100" }
      else
      {

        #print "$taxa[$taxon] $taxa[$taxon2]\n";
        if($ANIb_matrix{$taxa[$taxon]}{$taxa[$taxon2]})
        {

          #print scalar(@{$ANIb_matrix{$taxa[$taxon]}{$taxa[$taxon2]}})."\n";
          printf(ANIBMATRIX "\t%1.2f",
            calc_mean($ANIb_matrix{$taxa[$taxon]}{$taxa[$taxon2]}))
        }
        elsif($ANIb_matrix{$taxa[$taxon2]}{$taxa[$taxon]})
        {

          #print scalar(@{$ANIb_matrix{$taxa[$taxon2]}{$taxa[$taxon]}})."\n";
          printf(ANIBMATRIX "\t%1.2f",
            calc_mean($ANIb_matrix{$taxa[$taxon2]}{$taxa[$taxon]}))
        }
        else{ print ANIBMATRIX "\tNA" }
      }
    } print ANIBMATRIX "\n";
  }

  close(ANIBMATRIX);

  print "\n# average_identity_matrix_file = ".short_path($ANIb_matrix_file,$pwd)."\n";
  if($do_features){ print "# NOTE: matrix computed on blastn results\n" }
  else{ print "# NOTE: matrix computed on blastp results\n" }
}

# 5) produce clusters of intergenic sequences if requested
if($do_intergenic && !$total_clustersOK)
{
  if($total_genbankOK != $n_of_taxa)
  {
    die "# EXIT: cannot compile intergenic clusters as not all input GenBank files are valid ($total_genbankOK != $n_of_taxa)\n";
  }

  my ($ref_left_neighb,$ref_right_neighb,$ref_taxon,$n_of_clusterOK,$cluster_FNA,%intergenes);
  my (@left_orth_cluster,@right_orth_cluster,@ref_cluster,@other_cluster,$interg,$ref_interg);
  my ($ref_left_strand,$ref_right_strand,$ref_strand,$other_strand,$ref_interg_seq,$str1,$str2);
  my (@left_id_cluster,@right_id_cluster,@ref_id_cluster);
  my ($n_of_interg_clusters,$left_cluster_name,$right_cluster_name,$cluster_name) = (0);
  my ($interg_cluster_list_file,$interg_cluster_file,$cluster_FNA_counter,$interg_FNAresultsDIR);

  print "\n# looking for valid clusters of intergenic regions (n_of_taxa=$min_cluster_size)...\n";
  print "# parameters: MININTERGENESIZE=$MININTERGENESIZE MAXINTERGENESIZE=$MAXINTERGENESIZE INTERGENEFLANKORF=$INTERGENEFLANKORF\n\n";

  $interg_FNAresultsDIR     = $newDIR ."/".$output_mask."intergenic$MININTERGENESIZE\_$MAXINTERGENESIZE\_$INTERGENEFLANKORF\_";
  $interg_cluster_list_file = $interg_FNAresultsDIR . ".cluster_list";

  open(INTERGCLUSTERLIST,">$interg_cluster_list_file") || die "# EXIT : cannot create $interg_cluster_list_file\n";

  # 5.0) read FNA intergene files (this step requires some RAM)
  foreach $taxon (@taxa)
  {
    $intergenes{$taxon} = read_FASTA_file_array($intergenic_FNA_files{$taxon});
  }

  # 5.1) loop through intergenes of reference genome
  my $reference_intergenes = $intergenes{$taxa[$reference_proteome]};
  foreach $seq (0 .. $#{$reference_intergenes})
  {
    #>intergenic242|NC_005956|coords:1917497..1917790|length:294|neighbours:GI:49476308(-1),GI:49476309(1)|[Bartonella ...]
    # this is reference intergene
    if($reference_intergenes->[$seq][NAME] =~ m/neighbours:(.+?)\((\S+)\),(.+?)\((\S+)\)/)
    {
      @left_orth_cluster=@right_orth_cluster=@left_id_cluster=@right_id_cluster=();
      $left_cluster_name = $right_cluster_name = '';
      ($ref_left_neighb,$ref_left_strand,$ref_right_neighb,$ref_right_strand) = ($1,$2,$3,$4);

      # 5.1.1) check whether both intergene neighbors (left&right genes) are included
      # in orthologus protein clusters (OPCs), calculated in step #4
      foreach $cluster (keys(%GIclusters))
      {
        if(grep(/$ref_left_neighb/,@{$GIclusters{$cluster}}))
        {
          @left_orth_cluster = @{$GIclusters{$cluster}};
          @left_id_cluster = @{$idclusters{$cluster}};
          $left_cluster_name = $cluster;
        }
        if(grep(/$ref_right_neighb/,@{$GIclusters{$cluster}}))
        {
          @right_orth_cluster = @{$GIclusters{$cluster}};
          @right_id_cluster = @{$idclusters{$cluster}};
          $right_cluster_name = $cluster;
        }
        last if(@left_orth_cluster && @right_orth_cluster);
      }
      next if(!@left_orth_cluster || !@right_orth_cluster);

      # 5.1.2) both OPCs might have different size, choose smallest as @ref_cluster
      if(scalar(@left_orth_cluster) <= scalar(@right_orth_cluster))
      {
        @ref_cluster = @left_orth_cluster;
        @ref_id_cluster = @left_id_cluster;
        @other_cluster = @right_orth_cluster;
        $ref_strand = $ref_left_strand;
        $other_strand = $ref_right_strand;
        $cluster_name = "$left_cluster_name-$right_cluster_name";
      }
      else
      {
        @ref_cluster = @right_orth_cluster;
        @ref_id_cluster = @right_id_cluster;
        @other_cluster = @left_orth_cluster;
        $ref_strand = $ref_right_strand;
        $other_strand = $ref_left_strand;
        $cluster_name = "$right_cluster_name-$left_cluster_name";
      } #print "## $ref_strand $other_strand $cluster_name\n@ref_cluster\n@other_cluster)\n\n";

      # 5.1.3) check whether remaining OPCs members are as well neighbors of an intergene
      ($n_of_clusterOK,$cluster_FNA,$cluster_FNA_counter) = (0,'',1);
      foreach $orth (@ref_cluster)
      {

        # 5.1.3.1) find out intergenes flanked by this $orth gene (up to 2)
        my (@ref,@ref_seq);
        $ref_interg = $ref_interg_seq = '';
        TAXA: foreach $taxon (@taxa)
        {
          foreach $interg (@{$intergenes{$taxon}})
          {
            if($interg->[NAME] =~ /$orth/)
            {
              $ref_interg = $interg->[NAME];
              $ref_interg_seq = $interg->[SEQ];
              push(@ref,$ref_interg);
              push(@ref_seq,$ref_interg_seq);
              if(scalar(@ref)>1){ last TAXA }
            }
          }
        }
        next if(!@ref); # this should never happen

        # 5.1.3.2) find out whether any $other_orth gene in @other_cluster
        # is an intergenic neighbor of $orth
        OTHER: foreach my $other_orth (@other_cluster)
        {
          for(my $i=0;$i<scalar(@ref);$i++)
          {
            $ref_interg = $ref[$i];
            $ref_interg_seq = $ref_seq[$i];
            if($ref_interg =~ /$other_orth/)
            {

              # check that neighbors $orth,$other_orth and ref intergene
              # have equal strand orientation
              if($ref_interg =~ /$orth\((\S+)\),$other_orth\((\S+)\)/)
              {
                ($str1,$str2) = ($1,$2);
                if($str1 eq $ref_strand && $str2 eq $other_strand)
                {
                  $cluster_FNA .= ">$cluster_FNA_counter | $ref_interg\n";
                  $cluster_FNA .= "$ref_interg_seq\n";
                  $n_of_clusterOK++;
                  $cluster_FNA_counter++;

                  #print ">1 $orth $other_orth $str1 $str2 $ref_interg\n";
                }
              }
              elsif($ref_interg =~ /$other_orth\((\S+)\),$orth\((\S+)\)/)
              {
                ($str2,$str1) = ($1,$2);
                if($str1 ne $ref_strand && $str2 ne $other_strand) # strand can only be 1 or -1
                {
                  $ref_interg =~ s/\n/ (reverse complement)\n/;
                  $cluster_FNA .= ">$cluster_FNA_counter | $ref_interg\n";
                  $ref_interg_seq =~ tr/acgtACGT/tgcaTGCA/;
                  $ref_interg_seq = reverse($ref_interg_seq);
                  $cluster_FNA .= "$ref_interg_seq\n";
                  $n_of_clusterOK++; #print ">2 $other_orth $orth $str2 $str1 $ref_interg\n";
                  $cluster_FNA_counter++;
                }
              }
              last OTHER;
            }
          }
        }
      }

      # 5.1.4) produce clusters of orthologous intergenes
      if($n_of_clusterOK == scalar(@ref_cluster))
      {

        # check taxa present in this cluster
        my @interg_taxon_names;
        foreach my $id (@ref_id_cluster)
        {
          push(@interg_taxon_names,$gindex2[$id])
            if(!grep(/$gindex2[$id]/,@interg_taxon_names));
        }
        next if(scalar(@interg_taxon_names) < $min_cluster_size);

        $n_of_interg_clusters++;
        $cluster_name = $n_of_interg_clusters.'-'.$cluster_name;

        # create output directory if necessary
        if(!-e $interg_FNAresultsDIR){ mkdir($interg_FNAresultsDIR); }

        $interg_cluster_file = $interg_FNAresultsDIR ."/$cluster_name.fna";
        open(INTERGCLUSTER,">$interg_cluster_file") || die "# EXIT : cannot create $interg_cluster_file\n";
        print INTERGCLUSTER $cluster_FNA;
        close(INTERGCLUSTER);

        # print cluster info
        $annot = "cluster $cluster_name size=$n_of_clusterOK taxa=".
          scalar(@interg_taxon_names)." dnafile: $cluster_name.fna\n";
        print $annot if($PRINTCLUSTERSSCREEN);

        print INTERGCLUSTERLIST $annot;
        foreach $taxon (@interg_taxon_names){ print INTERGCLUSTERLIST ": $taxon\n"; }
      }
    }
  }

  if(!$n_of_interg_clusters){ print "\n# number_of_intergenic_clusters = $n_of_interg_clusters\n"; }
  else
  {
    print "\n# number_of_intergenic_clusters = $n_of_interg_clusters\n# intergenic_cluster_list = ".
      short_path($interg_cluster_list_file,$pwd)."\n" .
      "# intergenic_cluster_directory = ".short_path($interg_FNAresultsDIR,$pwd)."\n";
  }

  close(INTERGCLUSTERLIST);
}

# 6) clean and close stuff
if($saveRAM)
{
  untie(@sequence_data);untie(@sequence_prot);untie(@sequence_dna);
  unlink($sequence_data_filename,$sequence_prot_filename,$sequence_dna_filename); # not guaranteed

  if($do_PFAM)
  {
    untie (%pfam_hash);
    unlink ($pfam_data_filename);
  }
}

# 7) exit and show time
my $end_time = new Benchmark();
print "\n# runtime: ".timestr(timediff($end_time,$start_time),'all')."\n";
print "# RAM use: ".calc_memory_footprint()."\n";
exit(0);

#################################################################################

# based on http://www.unix.org.ua/orelly/perl/cookbook/ch04_18.htm
# generates a random permutation of @array in place
sub fisher_yates_shuffle
{
  my $array = shift;
  my ($i,$j,$array_string);

  for ($i = @$array; --$i; )
  {
    $j = int(rand($i+1));
    next if $i == $j;
    @$array[$i,$j] = @$array[$j,$i];
  }

  return join('',@$array);
}

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

# Split clusters of orths generated by OMCL, COG or PRND according to Pfam domain composition
# Updates $ref_orth_taxa and returns \%orthologues
# Uses globals: %pfam_hash, @gindex2
sub split_Pfam_clusters
{
  my ($ref_hash_orths, $ref_orth_taxa, $verbose) = @_;

  my ($pfam,$orth,$subcluster,%orthologues);
  my ($n_of_original_clusters,$n_of_new_clusters) = (0,0);

  print "# splitting clusters by Pfam domain composition\n";

  foreach $cluster (keys(%$ref_hash_orths))
  {

    # check number of Pfam domain strings in this cluster
    my (%Pfam_strings);

    $pfam = "$pfam_hash{$cluster}" || '';
    push(@{$Pfam_strings{$pfam}},$cluster);
    foreach $orth (@{$ref_hash_orths->{$cluster}})
    {
      $pfam = "$pfam_hash{$orth}" || '';
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

  printf("# created %d new clusters\n",$n_of_new_clusters-$n_of_original_clusters);

  return \%orthologues;
}

sub short_path
{
  my ($fullpath,$base) = @_;
  $fullpath =~ s/\/\//\//g;
  $fullpath =~ s/$base//;
  return $fullpath;
}

sub factorial
{
  my $max = int($_[0]);
  my $f = 1;
  for (2..$max) { $f *= $_ }
  return $f;
}

sub calc_mean
{
  my ($ref_args) = @_;
  my $mean = 0;
  foreach (@$ref_args) { $mean += $_ }
  return $mean / scalar(@$ref_args);
}

sub calc_std_deviation
{
  my ($ref_args) = @_;
  my $mean = calc_mean($ref_args);
  my @squares;

  foreach (@$ref_args){ push (@squares, ($_ **2)) }
  return sqrt( calc_mean(\@squares) - ($mean ** 2));
}

# works only in Linux systems, which have 'ps'
sub calc_memory_footprint
{
  my ($pid,$KB) = ($$,0);
  my $ps = `ps -o rss $pid 2>&1`; #`ps -o rss,vsz $pid 2>&1`
  while($ps =~ /(\d+)/g){ $KB += $1 }

  if($KB){ return sprintf("%3.1f MB",$KB/1024); }
  else
  {
    return "WARNING: cannot read memory footprint:\n$ps";
  }
}

# uses globasl: $0
sub warn_missing_soft
{
  my ($soft) = @_;
  print "# cannot run $soft as required software is not in place,\n";
  die "# EXIT : run $0 -v to check all required binares and data sources\n";
}

# checks whether SGE binaries are present and available in the system
# uses global $SGEPATH
sub cluster_is_available
{
  my $path = $SGEPATH;
  if($path eq '')
  {
    $path = `which qsub`;
    chomp($path);
    if($path eq '' || $path =~ /no qsub in/)
    {
      return 0;
    }
    else
    {
      $path =~ s/\/qsub$//;
      $SGEPATH = $path;
    }
  }

  for my $bin ('qstat','qsub','qdel')
  {
    my $output = `$path/$bin -help 2>&1`;
    if(!$output || $output !~ /usage:/){ return 0 }
  }

  return 1;
}

#sub estimate_RAM{
#  my ($n_of_sequences) = @_;
#  my ($use,$physicalRAM) = (0,0);
#  if(-s "/proc/meminfo") # will work only in linux systems
#  {
#    open(MEMINFO,"/proc/meminfo") || warn "# estimate_RAM : cannot read /proc/meminfo\n";
#    while(<MEMINFO>){ if(/MemTotal:\s+(\d+)/){ $physicalRAM = sprintf("%1.0f",$1/1024) } }
#    close(MEMINFO);
#  }
#  # bencFeb2012
#  #$use = sprintf("%1.0f",(0.030 * $n_of_sequences) + 130);
#  # bench Jun2013
#  $use = sprintf("%1.0f",(0.060 * $n_of_sequences) - 1700);
#  #if($isMCL){ $use = sprintf("%1.0f",(0.030 * $n_of_sequences) + 192) }
#  #else{ $use = sprintf("%1.0f",(0.027 * $n_of_sequences) + 122) }
#  if($physicalRAM && $use > $physicalRAM)
#  {
#    $use = "$use Mb (larger than physical RAM, consider using -s option)";
#  }
#  elsif($use > 0){ $use = "$use Mb" }
#  else{ $use = 0 }
#  return $use;
#}

# submits jobs to a SGE queue, this subroutine should be edited for other cluster systems
sub submit_cluster_job
{
  my ($jobname,$command,$outfile,$dir,$pathbin,$queuesets,$queuetime,$ref_cluster_PIDs) = @_;

  my $qsubcommand = " -N n$jobname -j y -o $outfile -S /bin/bash";
  my $qPID = `$pathbin/qsub $queuesets $qsubcommand <<EOF
        cd $dir
        $command
EOF`;

  if($qPID =~ /^(\d+)\./ || $qPID =~ /^Your job (\d+)/){ $qPID = $1 }
  $ref_cluster_PIDs->{$qPID}{'command'} = $qsubcommand;
  $ref_cluster_PIDs->{$qPID}{'executable'} = $command;
  $ref_cluster_PIDs->{$qPID}{'status'} = 'sent';
  sleep($queuetime);
}

# waits until launched cluster jobs terminate and handles failed jobs
sub check_cluster_jobs
{
  my ($dir,$pathbin,$queuetime,$waittime,$ref_PIDs) = @_;

  my ($waiting,$qPID,$newqPID,$qout) = (1);

  while($waiting)
  {
    $waiting=0;
    foreach $qPID (sort {$a<=>$b} (keys(%$ref_PIDs)))
    {
      next if($ref_PIDs->{$qPID}{'status'} eq 'deleted');
      $qout = `$pathbin/qstat | grep $qPID`;
      if($qout)
      {
        if($qout =~ /\s+$SGEERRORSTATUS\s+/)
        {

          # resubmit failed jobs
          $newqPID = `$pathbin/qsub $ref_PIDs->{$qPID}{'command'} <<EOF
            cd $dir
            $ref_PIDs->{$qPID}{'executable'}
EOF`;

          #Your job 108381 ("n_Bartonella_quintana_Toulouse.gbk.fasta") has been submitted
          #edit this regular expression if needed to suit your SGE cluster
          if($newqPID =~ /^(\d+)\./ || $newqPID =~ /^Your job (\d+)/){ $newqPID = $1 }
          $ref_PIDs->{$newqPID}{'command'} = $ref_PIDs->{$qPID}{'command'};
          $ref_PIDs->{$newqPID}{'executable'} = $ref_PIDs->{$qPID}{'executable'};
          $ref_PIDs->{$newqPID}{'status'} = 'sent';
          sleep($queuetime);

          # remove failed job
          system("$pathbin/qdel $qPID");
          $ref_PIDs->{$qPID}{'status'} = 'deleted';
          print "# check_cluster_jobs: deleted job $qPID , resubmitted as $newqPID\n";
          $waiting++;
        }
        else{ $waiting++; last; }
      }
    }
    if($waiting)
    {
      print "# check_cluster_jobs: waiting ...\n";
      sleep $waittime;
    }
  }
}
