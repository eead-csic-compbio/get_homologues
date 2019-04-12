#!/usr/bin/env perl

# 2016-8 Carlos P Cantalapiedra (1), Bruno Contreras-Moreira (1) and Pablo Vinuesa (2):
# 1: http://www.eead.csic.es/compbio (Laboratory of Computational Biology, EEAD/CSIC, Spain)
# 2: http://www.ccg.unam.mx/~vinuesa (Center for Genomic Sciences, UNAM, Mexico)

# This script checks whether some input clusters, produced by get_homologues.pl or get_homologues-est.pl,
# are enriched in Pfam domains when compared to previously Pfam-annotated sequences

# External dependencies: R , calls function fisher.test

$|=1;

use strict;
use warnings;
use Getopt::Std;
use File::Temp qw/tempfile/;
use List::Util 'sum';
use FindBin '$Bin';
use lib "$Bin/lib";
use lib "$Bin/lib/bioperl-1.5.2_102/";
use phyTools;
use marfil_homology;

my $VERBOSE = 0;  # set to 1 to see verbose messages
my @FEATURES2CHECK = ('EXE_R'); # cannot check EXE_PARS this way as it stalls
my @valid_tests = qw( greater two.sided less );
my $ADJUST = 'fdr'; # 'bonferroni' also supported

my (%opts,$INP_dir,$INP_clusterdir,$INP_clusterlist,$reference_proteome_string);
my ($INP_prot,$INP_test,$INP_pcutoff,$INP_est,$INP_seqs,$INP_fasta) = (1,'greater',0.05,0,0,0);

my ($id,$full_id,$cluster,$taxon,$pvalue,$adj,$pfam,$ref_pfam1,$ref_pfam2,$ref_full_text);
my (@clusters,%exp_clusters,%cluster_full_ids,@control_ids,@experiment_ids);
my (%cluster,%id2cluster,%fullid2id,%fasta,%fasta_annot);

getopts('hnsed:c:x:r:t:p:f:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-d directory with get_homologues[-est] output     (required, example: /path/to/data_homologues[-est])\n";
  print "-c path to directory with FASTA cluster files     (required, taken as 'control' set, sitting probably inside -d)\n";
  print "-x file with list of 'experiment' clusters        (required, should be a subset of files in -c)\n";
  print "-n use nucleotide sequence .fna clusters          (optional, uses .faa protein sequences by default)\n";
  print "-r reference .gbk/.faa/.fna file or \"taxon name\"  (optional, by default gets Pfam frequencies from all taxons)\n";
  print "-t type of test                                   (optional, default: -t greater, [".join(',',@valid_tests)."]\n";
  print "-p $ADJUST adjusted p-value threshold                 (optional, default: -p $INP_pcutoff\n";
  print "-e cluster sequences are ESTs                     (optional, use with pre-computed Pfam mappings of deduced CDSs)\n";
  print "-f make FASTA file with experiment sequences      (optional, example: -f experiment.fasta)\n";
  print "-s take individual sequences as annotation units  (optional, by default Pfam domains are counted once per cluster)\n\n"; 
  exit(0); 
}

if(defined($opts{'d'}))
{  
  $INP_dir = $opts{'d'}; 
}
else{ die "# EXIT : need a -d directory\n"; }

if($opts{'n'}){ $INP_prot = 0 }

if(defined($opts{'c'}) && defined($opts{'x'}))
{  
  # read cluster names from folder
  $INP_clusterdir = $opts{'c'}; 
  opendir(CLUSTERDIR,$INP_clusterdir) || die "# $0 : cannot list $INP_clusterdir\n";
  
  my @tmpclusters = readdir(CLUSTERDIR);
  
  if($INP_prot)
  {
    @clusters = grep{/\.faa$/} grep{!/pangenome_matrix/} @tmpclusters;
  }
  else
  {
    @clusters = grep{/\.fna$/} grep{!/pangenome_matrix/} @tmpclusters;
  }  
 
  closedir(CLUSTERDIR);
  
  die "# $0 : cannot find valid clusters in $INP_clusterdir,\n".
    "# make sure you set -n with nucleotide clusters\n" if(!@clusters);
   
  # read list of experiment clusters
  $INP_clusterlist = $opts{'x'}; 
  open(CLUSTERLIST,$INP_clusterlist) || die "# $0 : cannot read $INP_clusterlist\n";
  while(<CLUSTERLIST>)
  {
    next if(/^#/);
    $exp_clusters{(split)[0]} = 1; #9006_TR29928-c0_g1_i1.fna
  }
  close(CLUSTERLIST);
  
  die "# $0 : cannot parse 'experiment' clusters in $INP_clusterlist,\n".
    "# make sure you set -n with nucleotide clusters\n" if(!keys(%exp_clusters));
}
else{ die "# EXIT : need -c path and -x file\n"; }

if($opts{'r'}){ $reference_proteome_string = $opts{'r'}; }
else{ $reference_proteome_string = '' }

if($opts{'t'})
{ 
  $INP_test = $opts{'t'};
  if(!grep(/$INP_test/,@valid_tests))
  {
    die "# $0 : valid test types are: ".join(',',@valid_tests)."\n";
  }  
}

if($opts{'p'} && $opts{'p'} > 0.0 && $opts{'p'} <= 1.0)
{ 
  $INP_pcutoff = $opts{'p'};
}

if($opts{'e'}){ $INP_est = 1 }

if($opts{'s'}){ $INP_seqs = 1 }

if($opts{'f'}){ $INP_fasta = $opts{'f'} }

printf("# %s -d %s -c %s -x %s -n %d -s %d -e %d -t %s -p %1.2f -r %s -f %s\n\n",
  $0,$INP_dir,$INP_clusterdir,$INP_clusterlist,!$INP_prot,$INP_seqs,$INP_est,
  $INP_test,$INP_pcutoff,$reference_proteome_string,$INP_fasta);

##########################################################################

# 0) check that R can be called
check_installed_features(@FEATURES2CHECK);
if(!feature_is_installed('R'))
{
  print "\n# WARNING : this script requires the software R, available from http://www.r-project.org\n";
  exit;
}


## 1) check whether precomputed Pfam-assignments are in place
constructDirectory($INP_dir); # this sets $pfam_file and $p2ofilename
if(!-s $pfam_file)
{
  die "# ERROR : cannot find $pfam_file (generated by get_homologues[-est].pl -D), exit\n";
}


## 2) parse full sequence ids from input clusters 
print "# parsing clusters...\n";

foreach my $cluster (@clusters)
{
  if($INP_est) # Takes first non-blank tag + [taxon] in FASTA header to identify sequences
  {
    open(CLUSTER,"$INP_clusterdir/$cluster") || 
      die "# $0 : cannot read list $INP_clusterdir/$cluster\n"; 
    while(<CLUSTER>)
    {
      chomp;
      if(/^>(\S+).*?\[(\S+?)\.fna/)
      { 
        #>TR15703|c0_g1_i1 evidence:transdecoder   [Alexis.trinity.fna_minl50_eval1e-05.cds.fna.bz2] | aligned:13-498 (498)
        #>TR32976|c0_g1_i2 len=3138 path=[6380:0-378...] [Alexis.trinity.fna.bz2]
        #evidence:transdecoder.blastx match:sp|Q65XV7|RFA1C_ORYSJ
        ($full_id,$taxon) = ($1,$2); 
        
        if($reference_proteome_string && $taxon !~ /$reference_proteome_string/)
        {
          $full_id = '';
          next;
        }  
        
        $full_id = $full_id." [$taxon]"; #TR30768|c0_g1_i1SBCC073_fLF.Trinity
      
        if($cluster_full_ids{$full_id})
        {
          warn "# WARNING: skipping sequence: $full_id , seems to be duplicated:\n$cluster{$full_id},$cluster\n";
        }
        else
        { 
          $cluster_full_ids{$full_id} = $_; 
          $cluster{$full_id} = $cluster; 
          
          if($INP_fasta && /(evidence:\S+ match:\S+)/ || /(evidence:\S+)/)
          {
            $fasta_annot{$full_id} = $1;
          }
        }  
      } 
      else
      {
        if($INP_fasta && $full_id ne '')
        {
          $fasta{$full_id} .= $_; 
        }
      } 
    }
    close(CLUSTER);
  }
  else
  {
    open(CLUSTER,"$INP_clusterdir/$cluster") || 
      die "# $0 : cannot read list $INP_clusterdir/$cluster\n";
    while(<CLUSTER>)
    {
      chomp;
      if(/^>(.*)?\|neighbours/ || /^>(.*)? \| aligned/ || /^>(.*)/)
      { 
        #>gi|116515001|ref|YP_802630.1| DapE [Buchnera aphidicola str. Cc (Cinara cedri)] | aligned:1-375 (376)
        
        if($reference_proteome_string && $_ !~ /$reference_proteome_string/)
        {
          $full_id = '';
          next;
        }
        
        $full_id = $1;
        
        if($cluster_full_ids{$full_id})
        {
          warn "# WARNING: skipping sequence $full_id , seems to be duplicated:\n$cluster{$full_id},$cluster\n";
        }
        else
        { 
          $cluster_full_ids{$full_id} = $_;
          $cluster{$full_id} = $cluster; 
        }   
      }
      else
      {
        if($INP_fasta && $full_id ne '')
        {
          $fasta{$full_id} .= $_; 
        }
      } 
    }
    close(CLUSTER);
  }
}

if(scalar(keys(%cluster_full_ids)) == 0)
{
  printf("# ERROR : did not extract any sequences from %d clusters, exit\n",
    scalar(@clusters));
  exit;
}

printf("# %d sequences extracted from %d clusters\n\n",
  scalar(keys(%cluster_full_ids)),scalar(@clusters));
  


## 3) get short ids for selected sequences
open(P2O,$p2ofilename) || 
  die "# ERROR : cannot find $p2ofilename (generated by get_homologues[-est]), exit\n";
while(<P2O>)
{
  #1,Buch_aph_APS.faa,gi|10957100|ref|NP_057962.1| anthranilate synthase component I [Buchnera aphidicola str. APS (Acyrthosiphon pisum)]
  #3073,Alexis.trinity.fna_minl50_eval1e-05.cds.fna.bz2,TR15399|c0_g1_i1 evidence:transdecoder  
  
  # skip non-reference sequences if requested
  next if($reference_proteome_string && $_ !~ m/$reference_proteome_string/);
  
  chomp;
  ($id,$taxon,$full_id) = split(/,/,$_,3);
  
  if($INP_est)
  {
    if($taxon =~ m/(\S+?)\.fna/){ $taxon = $1 } 
    $full_id = (split(/\s+/,$full_id))[0]." [$taxon]"; #print "$full_id\n";
  }   
  else
  {
    if($full_id =~ m/(.*)?\|neighbours/){ $full_id = $1 }
    elsif($full_id =~ m/(.*)?\| aligned/){ $full_id = $1 }
  }
 
  # find out cluster to which this id belongs
  if($cluster{$full_id})
  {
    $cluster = $cluster{$full_id};
  }
  else{ next } # skip unclustered sequences to avoid biass
  
  # store 'experiment' sequence ids
  if($exp_clusters{$cluster})
  {
    push(@experiment_ids,$id);
    
    # save equivalences between ids of experiment sequences
    if($INP_fasta)
    {
      $fullid2id{$full_id} = $id; 
    }  
  }
  
  # store 'control' ids 
  push(@control_ids,$id); # naturally sorted
  
  # save table of ids and clusters
  $id2cluster{$id} = $cluster;
}
close(P2O);

if(scalar(@experiment_ids) == 0)
{
  print "# ERROR : cannot match ids of 'experiment' sequences, exit\n";
  exit;
}

printf("# total experiment sequence ids = %d\n",scalar(@experiment_ids));
printf("# total control    sequence ids = %d\n\n",scalar(@control_ids)); 

foreach $full_id (keys(%cluster_full_ids))
{
  print "# cannot find id of: $full_id\n" if(!defined($cluster_full_ids{$full_id}) && $VERBOSE);
}


## 4) get Pfam annotations of selected sequences (experiment) and control
if($INP_seqs)
{
  ($ref_pfam1,$ref_pfam2,$ref_full_text) = parse_Pfam_freqs($pfam_file,\@experiment_ids,\@control_ids); 
}
else
{
  ($ref_pfam1,$ref_pfam2,$ref_full_text) = parse_Pfam_freqs($pfam_file,\@experiment_ids,\@control_ids,\%id2cluster);
}  

if(scalar(keys(%$ref_pfam1)) != scalar(keys(%$ref_pfam2)))
{
  printf("# ERROR : experiment (%d) and control (%d) Pfam annotation tables differ in length, exit\n",
    scalar(keys(%$ref_pfam1)),scalar(keys(%$ref_pfam2)));
  exit;
}

if($INP_fasta)
{
  # fill $marfil_homology::pfam_hash
  construct_Pfam_hash($pfam_file);
  
  # print experiment sequences cluster by cluster
  my ($mean_length,$seqs_per_cluster,$total,%clusters) = (0,0,0);
  open(OUTFASTA,">",$INP_fasta) || die "# ERROR : cannot create file $INP_fasta\n";
  foreach $full_id (sort {$cluster{$a} cmp $cluster{$b}} keys(%fasta))
  { 
    next if(!$exp_clusters{$cluster{$full_id}});
    
    printf( OUTFASTA ">%s cluster=%s %s %s\n%s\n",
      $full_id,$cluster{$full_id},
      $fasta_annot{$full_id} || '', 
      $pfam_hash{$fullid2id{$full_id}},
      $fasta{$full_id});
      
    $mean_length += length($fasta{$full_id});
    $clusters{$cluster{$full_id}}++;
    $total++;
  } 
  close(OUTFASTA);
  
  $mean_length = sprintf("%1.1f",$mean_length/$total);
  $seqs_per_cluster = sprintf("%1.2f",$total/scalar(keys(%clusters)));
  
  print "\n# created FASTA file: $INP_fasta\n\n";
  print "# sequences=$total mean length=$mean_length , seqs/cluster=$seqs_per_cluster\n\n"; 
}

my ($fh1,$filename1) = tempfile(SUFFIX => '.dat', UNLINK => 1); # experiment (query in R)
my ($fh2,$filename2) = tempfile(SUFFIX => '.dat', UNLINK => 1); # control (reference in R)
my ($fh3,$filename3) = tempfile(SUFFIX => '.dat', UNLINK => 1); # enrichment results
my ($fh4,$filename4) = tempfile(SUFFIX => '.dat', UNLINK => 1); # sort output

foreach $pfam (sort (keys(%$ref_pfam1)))
{
  print $fh1 "$pfam\t$ref_pfam1->{$pfam}\n";
}

foreach $pfam (sort (keys(%$ref_pfam2)))
{
  print $fh2 "$pfam\t$ref_pfam2->{$pfam}\n";
}

close($fh3);
close($fh2);
close($fh1);

## 5) calculate enrichment
my $Rparams = '';
if(!$VERBOSE){ $Rparams = '-q 2>&1 > /dev/null' }

open(RSHELL,"|R --no-save $Rparams ") || die "# $0 : cannot call R: $!\n";
print RSHELL<<EOF;

# globals
direction="$INP_test"
multitest="$ADJUST"
verbose=F

# parse data
que_data=read.csv(file="$filename1", sep="\t", header=FALSE);
que_rows=nrow(que_data);
ref_data=read.csv(file="$filename2", sep="\t", header=FALSE);
ref_rows=nrow(ref_data); 

# uses globals que_total, ref_total
test <- function(x){
  que_id=x[1];
  que_value=as.numeric(x[2]);
  ref_value=as.numeric(ref_data[ref_data\$V1==que_id,2]);

  if (length(ref_value)==0){
    ref_value=0;
  }

  if(verbose==T){
    cat(paste(que_id, "\n"), file=stderr());
  }

  values=c(que_value, ref_value, que_total, ref_total);
  input_matrix=matrix(values, nrow = 2,
    dimnames=list(c("exp", "control"), c("Pfam", "total")));

  if(verbose==T){
    cat(paste(input_matrix, "\n"), file=stderr());
  }

  fisher_htest=fisher.test(input_matrix, alternative=direction);

  if(verbose==T){
    cat(paste(fisher_htest, "\n"), file=stderr());
  }

  ret_value=c( que_id, fisher_htest\$p.value );

  return(ret_value);
}

print_pvalues <- function(x) {
  cat(paste(x[1],"\t",x[2],"\t",x[3],"\n"),file="$filename3",append=T);
}

que_total=sum(que_data[,2]);
ref_total=sum(ref_data[,2]);

pvalues=apply(que_data, 1, test);

## Multiple test adjustment
num_pvalues=as.numeric(pvalues[2,]);
adj_pvalues=p.adjust(num_pvalues, method=multitest);

result_pvalues=rbind(pvalues, adj_pvalues);
result_pvalues=t(result_pvalues); # rows to columns

apply(result_pvalues, 1, print_pvalues);
q()
EOF
close RSHELL;

## 6) print enrichment results

# count total assigned Pfam domains in both sets
my $total_exp = sum(values(%$ref_pfam1));
my $total_ctr = sum(values(%$ref_pfam2));

if($total_exp==0)
{
  die "# ERROR : zero parsed Pfam domains (experiment), please check -x list\n";
}

print "# fisher exact test type: '$INP_test'\n";
print "# multi-testing p-value adjustment: $ADJUST\n";
print "# adjusted p-value threshold: $INP_pcutoff\n\n";
printf("# total annotated domains: experiment=%d control=%d\n\n",$total_exp,$total_ctr);

system("$SORTBIN -n -k3,3 -k2,2 $filename3 > $filename4");
if($? != 0)
{
  die "# ERROR: failed running $SORTBIN -n -k3,3 -k2,2 $filename3 > $filename4 \n";
}

open(ENRICH,"<",$filename4) || die "# $0 : cannot read enrichment results file $filename4\n";
print "#PfamID\tcounts(exp)\tcounts(ctr)\tfreq(exp)\tfreq(ctr)\tp-value\tp-value(adj)\tdescription\n";
while(<ENRICH>)
{
  chomp;
  ($pfam,$pvalue,$adj) = split(/\s+/,$_);

  last if($adj > $INP_pcutoff);

  $ref_pfam1->{$pfam} |= 0;
  $ref_pfam2->{$pfam} |= 0;

  printf("%s\t%d\t%d\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%s\n",
    $pfam,
    $ref_pfam1->{$pfam},$ref_pfam2->{$pfam},
    $ref_pfam1->{$pfam}/$total_exp,$ref_pfam2->{$pfam}/$total_ctr,
    $pvalue,$adj,$ref_full_text->{$pfam});
}
close(ENRICH);

