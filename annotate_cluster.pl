#!/usr/bin/env perl

# 2016 Bruno Contreras-Moreira (1) and Pablo Vinuesa (2):
# 1: http://www.eead.csic.es/compbio (Laboratory of Computational Biology, EEAD/CSIC, Spain)
# 2: http://www.ccg.unam.mx/~vinuesa (Center for Genomic Sciences, UNAM, Mexico)

# This script produces a multiple alignment view of BLAST locally-aligned sequences clustered by get_homologues[-est].pl
# It does not necessarily conserved the original sequence order 

# Uses third-party software:
# MVIEW (https://github.com/desmid/mview) Brown NP, Leroy C, Sander C (1998) MView: A Web compatible database search or 
# multiple alignment viewer. Bioinformatics. 14 (4):380-381.

$|=1;

use strict;
use Getopt::Std;
use FileHandle;
use File::Temp qw/tempfile/;
use FindBin '$Bin';
use lib "$Bin/lib";
use lib "$Bin/lib/bioperl-1.5.2_102/";
use Bio::Seq;
use Bio::SeqUtils;
use phyTools;
use marfil_homology;

my @FEATURES2CHECK = ('EXE_BLASTP','EXE_BLASTN','EXE_FORMATDB','EXE_MVIEW','EXE_HMMPFAM');

my $DEFBLASTNTASK  = 'megablast';
my $DEFEVALUE      = 10; # default BLAST E-value
my $MINBLUNTBLOCK  = 100; # min alignment width with blunt ends
my $MAXSEQNAMELEN  = 30;

my ($INP_nucleotides,$INP_blunt,$do_PFAM,$INP_clusterfile,$INP_outfile,%opts) = (1,0,0,'','');

getopts('hDbPf:o:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-f input cluster FASTA file            (expects nucleotides)\n";        
  print   "-o output alignment file               (optional, produces FASTA format)\n";
  print   "-P sequences are peptides              (optional)\n";
  print   "-b blunt alignment borders             (optional, otherwise all bits aligned by BLAST are shown)\n";
  print   "-D match Pfam domains                  (optional, based on six-frame translation of longest seq)\n";
  exit(0);
}

print STDERR check_installed_features(@FEATURES2CHECK);

if(defined($opts{'f'})){  $INP_clusterfile = $opts{'f'}; }
else{ die "# EXIT : need -f cluster FASTA file\n"; }

if(defined($opts{'o'})){  $INP_outfile = $opts{'o'}; }

if(defined($opts{'P'})){  $INP_nucleotides = 0 }

if(defined($opts{'b'}))
{ 
	$INP_blunt = $MINBLUNTBLOCK;
}

if(defined($opts{'D'}))
{
  if(feature_is_installed('PFAM'))
  {
    $do_PFAM = 1;
  }
  else{ warn_missing_soft('PFAM') }
}

printf(STDERR "\n# %s -f %s -o %s -P %s -b %s -D %s\n",
  $0,$INP_clusterfile,$INP_outfile,$INP_nucleotides,$INP_blunt,$do_PFAM);

##########################################################################

my ($maxlength,$longest_seq,$length,$start,$end,$seq,$fh,$seqname) = (0,-1);
my ($command,$Pfam_string,$Pfam_full,@trash,%short_names);

## parse input cluster and get sequences
my $cluster_ref = read_FASTA_file_array($INP_clusterfile,$INP_nucleotides);
foreach $seq (0 .. $#{$cluster_ref})
{
  # shorten name
  if($cluster_ref->[$seq][NAME] =~ m/^(\S+).*?(\[\S+?\])$/ || 
    $cluster_ref->[$seq][NAME] =~ m/^(\S+).*?(\[\S+?\])\s\|/)
  {
    $seqname = "$1\_$2"; 
  }
  else
  {
    $seqname = (split(/\s+/,$cluster_ref->[$seq][NAME]))[0];
  }
  
  # shorten name, otherwise it might break mview downstream
  if(length($seqname) > $MAXSEQNAMELEN ){ $seqname = substr($seqname,0,$MAXSEQNAMELEN ) } 
  
  $short_names{$seq} = $seqname;

  # in case of missing sequences
  if(!$cluster_ref->[$seq][SEQ])
  { 
    print STDERR "# skip $seq : $cluster_ref->[$seq][NAME] due to missing sequence\n";
    next;
  }  

  # find longest aligned sequence to be used as reference | aligned:1-296 (296)
  if($cluster_ref->[$seq][NAME] =~ / \| aligned:(\d+)-(\d+) /)
  { 
    ($start,$end) = ($1,$2);
    $length = ($end-$start)+1;
  }
  else
  {
    $length = length($cluster_ref->[$seq][SEQ]);
  }
  
  # find longest [aligned] sequence
  if($maxlength == 0 || $length > $maxlength)
  {
    $maxlength = $length;
    $longest_seq = $seq;
  }
}

#printf(STDERR "\n# total   sequences: %d longest: %d\n",$#{$cluster_ref}+1,$longest_seq+1); 
printf(STDERR "\n# total   sequences: %d\n",$#{$cluster_ref}+1); 

if($#{$cluster_ref} == -1)
{
  warn "# Cannot read input sequences, exit. Please set -P if using a peptide cluster.\n";
  exit;
}

## Pfam-annotate longest sequence if required
if($do_PFAM)
{
  my ($fhpf,$filenamepf) = tempfile(UNLINK => 1); # Pfam input 
  my ($fhpo,$filenamepo) = tempfile(UNLINK => 1); # Pfam output
  my ($fhpr,$filenamepr) = tempfile(UNLINK => 1); # Pfam report

  if($INP_nucleotides)
  {
    # translate to 6 frames 
    # https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/
    my $seqobj = Bio::Seq->new (-seq=>$cluster_ref->[$longest_seq][SEQ],-alphabet=>'DNA');
    my $utils = new Bio::SeqUtils;
    my @seqs = $utils->translate_6frames($seqobj);
    foreach $seq (0 .. $#seqs)
    {
      print $fhpf ">$seq\n$seqs[$seq]->{'primary_seq'}{'seq'}\n";
    }
  }
  else
  {
    print $fhpf ">1\n$cluster_ref->[$longest_seq][SEQ]\n";
  }  
  
  close($fhpf);
  close($fhpo);
  close($fhpr);
  
  $command = format_HMMPFAM_command()." $filenamepf > $filenamepo ";
  system($command);
  if($? != 0)
  {
    die "# EXIT: failed while running Pfam search ($command)\n";
  }#else{ system("cat $filenamepo"); } # testing
  
  pfam_parse($filenamepr,($filenamepo));
  
  open(PFAMREPORT,$filenamepr);
  while(<PFAMREPORT>)
  {
    if(/^\d+\t(PF\S+)\t(.*?)\n/)
    {
      $Pfam_string .= $1;
      $Pfam_full   .= $2;
    }  
  }
  close(PFAMREPORT);
  
  printf(STDERR "\n# Pfam domains: %s\n",$Pfam_string); 
  printf(STDERR "# Pfam annotation: %s\n",$Pfam_full); 
}  
  
## align cluster sequences
# temp files, removed automatically on exit
my ($fhdb,$filenamedb) = tempfile(UNLINK => 1); # blast DB
my ($fhq,$filenameq)   = tempfile(UNLINK => 1); # blast query
my ($fhb,$filenameb)   = tempfile(UNLINK => 1); # blast output
my ($fha,$filenamea)   = tempfile(UNLINK => 1); # out alignment file

foreach $seq (0 .. $#{$cluster_ref})
{
  if($seq == $longest_seq){ $fh = $fhq }
  else{ $fh = $fhdb }
  
  print $fh ">$short_names{$seq}\n$cluster_ref->[$seq][SEQ]\n";
}

close($fha);
close($fhb);
close($fhq);
close($fhdb);

if($INP_nucleotides)
{
  executeFORMATDB($filenamedb,1,1); 
  $command = format_BLASTN_command_aligns($filenameq,$filenameb,$filenamedb,$DEFEVALUE,$DEFBLASTNTASK);
  push(@trash,$fhdb.'.nsq', $fhdb.'.nin', $fhdb.'.nhr');
}
else
{
  executeFORMATDB($filenamedb,0,1);
  $command = format_BLASTP_command_aligns($filenameq,$filenameb,$filenamedb,$DEFEVALUE);
  push(@trash,$fhdb.'.psq', $fhdb.'.pin', $fhdb.'.phr');
}

system($command);
if($? != 0)
{
  die "# EXIT: failed while running BLAST search ($command)\n";
}
else
{
  # clean BLAST tmp files
  unlink(@trash);
}

## format raw alignment  #-width 60 
$command = "$ENV{'EXE_MVIEW'} -in blast -out fasta $filenameb > $filenamea"; 
system($command); 

## extract SNPs and parsimony-informative sites, and trim alignment if required
my $align_ref = check_variants_FASTA_alignment($filenamea,!$INP_nucleotides,$INP_blunt);
  
printf(STDERR "# aligned sequences: %d width:   %d\n",$#{$align_ref}+1,length($align_ref->[0][SEQ]));  
  
## print final alignment 
if($INP_outfile)
{
  $fh = FileHandle->new(">$INP_outfile");
  if(!defined($fh))
  {
    die "# cannot create output file $INP_outfile\n";
  }
}
else{ $fh = *STDOUT }

foreach $seq (0 .. $#{$align_ref})
{
  print $fh ">$align_ref->[$seq][NAME]";
  
  if($do_PFAM){ print $fh " Pfam:$Pfam_full($Pfam_string)" }
  
  print $fh "\n$align_ref->[$seq][SEQ]\n";
} 

if($INP_outfile && $#{$align_ref})
{
  print STDERR "# alignment file: $INP_outfile\n"; 
}  



