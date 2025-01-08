#!/usr/bin/env perl

# 2016-25 Bruno Contreras-Moreira (1), Ruben Sancho (1,2) and Pablo Vinuesa (3):
# 1: httpa://www.eead.csic.es/compbio (Estacion Experimental Aula Dei-CSIC/Fundacion ARAID, Spain)
# 2: Grupo Bioflora, EPS, Universidad de Zaragoza, Spain
# 3: http://www.ccg.unam.mx/~vinuesa (Center for Genomic Sciences, UNAM, Mexico)

# This script produces a multiple alignment view of BLAST locally-aligned sequences clustered by
# get_homologues[-est].pl
# It works by aligning all sequences in the cluster to the longest/user-selected sequence.
# It does not necessarily conserve the original sequence order. 
# Any sequence stretches masked by BLAST will appear as Xs.

# Uses third-party software:
# MVIEW (https://github.com/desmid/mview) Brown NP, Leroy C, Sander C (1998) 
# MView: A Web compatible database search or multiple alignment viewer. 
# Bioinformatics. 14 (4):380-381.

$|=1;

use strict;
use warnings;
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

my $DEFBLASTNTASK = 'megablast';
my $DEFEVALUE     = 10; # default BLAST E-value
my $MINBLUNTBLOCK = 100; # min alignment width with blunt ends
my $MAXSEQNAMELEN = 60;
my $MAXMISMCOLLAP = 0; # natural, number of mismatches tolerated when collapsing 
my $MAXGAPSCOLLAP = 2; # natural, number of tolerated gaps in different places when collapsing 

my %CHARSFORMVIEW = ( '/' => '_fslash_' ); # chars that have to be removed before calling MVIEW
my %CHARSFROMMVIEW = ( '_fslash_' => '/' );

my ($INP_nucleotides,$INP_blunt,$do_PFAM,$INP_collapse) = (1,0,0,0);
my ($INP_clusterfile,$INP_outfile,$INP_ref_file,%opts) = ('','','');
my ($INP_unaligned,$INP_includeA,$INP_includeB) = (0,'','');

my $warning = <<'END_WARN';
WARNING1: Clusters of transcripts often contain a fraction of BLAST hits that do not match 
the longest sequence; instead, they align towards the 5' or 3' of other sequences and are 
not included in the produced cumulative multiple sequence alignment (MSA):
 
 -----------------            <= longest/reference sequence
              -------------
    -------------
 -----------
   ------------
                      ....    <= sequences not included in MSA
                       ..

WARNING2: The resulting multiple alignment is produced by MVIEW, which does not record
deletions in the longest/reference sequence. This means that an alignment like this:

 ------  -----------            <= longest/reference sequence
    ---..----------             <= .. fragment not included in MSA
 ------  -----
   ----  --------

will in fact be saved as:

-----------------            
   -------------             
-----------
  ------------

END_WARN

getopts('hDbPuf:o:r:c:A:B:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-f input cluster FASTA file          (expects nucleotides, aligns longest seq to rest of cluster)\n";        
  print   "-o output FASTA file                 (optional, by default prints alignment to standard output)\n";
  print   "-P sequences are peptides            (optional, uses BLASTP instead of BLASTN)\n";
  print   "-r reference sequence FASTA          (optional, aligns cluster sequences to this external seq)\n";
  print   "-D match Pfam domains                (optional, annotates longest seq, nucleotides on 6 frames)\n";
  print   "-u print unaligned sequences         (optional, flips seqs, handy for multiple alignments, skips option below)\n";
  print   "-b blunt alignment borders           (optional, also annotates SNPs/missense mutations and parsimony-informative sites)\n";
  print   "-A file with taxon names of group A  (optional, identifies private variants of group A vs 'rest')\n";
  print   "-B file with taxon names of group B  (optional, requires -A, group B is used as 'rest')\n";
  print   "-c collapse overlapping fragments    (optional, example -c 40 for overlaps >= 40 residues, requires -o,\n";
  print   "                                      this is useful to merge fragmented de-novo transcripts)\n\n";
  print   $warning;
  exit(0);
}

print STDERR check_installed_features(@FEATURES2CHECK);

if(defined($opts{'f'})){  $INP_clusterfile = $opts{'f'}; }
else{ die "# EXIT : need -f cluster FASTA file\n"; }

if(defined($opts{'o'}))
{  
  $INP_outfile = $opts{'o'}; 

  if(!defined($opts{'u'}) && defined($opts{'c'}) && $opts{'c'} > 0)
  {
    $INP_collapse = $opts{'c'};
  }
}

if(defined($opts{'r'})){  $INP_ref_file = $opts{'r'}; }

if(defined($opts{'D'}))
{
  if(feature_is_installed('PFAM'))
  {
    $do_PFAM = 1;
  }
  else{ warn_missing_soft('PFAM') }
}

if(defined($opts{'P'})){  $INP_nucleotides = 0 }

if(defined($opts{'u'})){  $INP_unaligned = 1 }
else
{
  if(defined($opts{'b'}))
  { 
    $INP_blunt = $MINBLUNTBLOCK;
  }

  if(defined($opts{'A'}))
  {
    $INP_includeA = $opts{'A'};

    if(defined($opts{'B'}))
    {
      $INP_includeB = $opts{'B'};
    }
  }

}

warn "\n# DEFBLASTNTASK=$DEFBLASTNTASK DEFEVALUE=$DEFEVALUE\n";
warn "# MINBLUNTBLOCK=$MINBLUNTBLOCK MAXSEQNAMELEN=$MAXSEQNAMELEN\n";
warn "# MAXMISMCOLLAP=$MAXMISMCOLLAP MAXGAPSCOLLAP=$MAXGAPSCOLLAP\n";

printf(STDERR "\n# %s -f %s -r %s -o %s -P %d -u %d -b %s -D %s -c %d -A %s -B %s\n",
  $0,$INP_clusterfile,$INP_ref_file,$INP_outfile,!$INP_nucleotides,
  $INP_unaligned,$INP_blunt,$do_PFAM,$INP_collapse,$INP_includeA,$INP_includeB);

##########################################################################

my ($maxlength,$longest_seq,$length,$start,$end,$seq,$fh,$seqname) = (0,-1);
my ($ext_seq,$ext_seqname,$taxon,$seqid,$align_seq,$char) = ('','');
my ($command,$Pfam_string,$Pfam_full) = ('','','');
my (@trash,%short_names,%intaxa,%id2taxon);
my (%included_input_filesA,%included_input_filesB,@listA,@listB);

## read reference external sequence if requested (take first seq only)
if($INP_ref_file ne "")
{
  my $external_ref = read_FASTA_file_array($INP_ref_file,$INP_nucleotides);
  foreach $seq (0 .. $#{$external_ref})
  {
    # shorten name
    $taxon = '';
    if($external_ref->[$seq][NAME] =~ m/^(\S+).*?\[(\S+?)\]$/ || 
      $external_ref->[$seq][NAME] =~ m/^(\S+).*?\[(\S+?)\]\s\|/ ||
      $external_ref->[$seq][NAME] =~ m/^(\S+).*?\[(.+?)\]\s\|/ ||
      $external_ref->[$seq][NAME] =~ m/^(\S+).*?\[(.+?)\]/)
    {
      $taxon = (split(/\.f/,$2))[0];
      $ext_seqname = "$1\[$taxon]";
      $ext_seqname =~ s/\s+/_/g;
      
      $intaxa{$taxon}++;
    }
    else{ $ext_seqname = (split(/\s+/,$external_ref->[$seq][NAME]))[0] }

    $ext_seq = $external_ref->[$seq][SEQ];
    last;
  }

  warn "# external sequence parsed: $ext_seqname\n";
}

## parse input cluster and get sequences
my $cluster_ref = read_FASTA_file_array($INP_clusterfile,$INP_nucleotides);

foreach $seq (0 .. $#{$cluster_ref})
{
  # shorten name
  $taxon = '';
  if($cluster_ref->[$seq][NAME] =~ m/^(\S+).*?\[(\S+?)\]$/ ||  
    $cluster_ref->[$seq][NAME] =~ m/^(\S+).*?\[(\S+?)\]\s\|/ ||
    $cluster_ref->[$seq][NAME] =~ m/^(\S+).*?\[(.+?)\]\s\|/ ||
    $cluster_ref->[$seq][NAME] =~ m/^(\S+).*?\[(.+?)\]/)
  {
    ($seqid,$taxon) = ($1,$2);

    $taxon = (split(/\.f/,$taxon))[0];
    $seqname = "$seqid\[$taxon]"; 
    $seqname =~ s/\s+/_/g;
    
    $intaxa{$taxon}++;
    
    # link id to taxon, as ids are less likely to be truncated in BLAST output
    $id2taxon{$seqid} = $taxon;
  }
  else
  {
    $seqname = (split(/\s+/,$cluster_ref->[$seq][NAME]))[0];
  }
 
  # shorten name and remove forbidden chars, otherwise it might break mview downstream
  if(length($seqname) > $MAXSEQNAMELEN )
  { 
    $seqname = substr($seqname,0,$MAXSEQNAMELEN); 
  } 
  
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

if(!defined($cluster_ref->[0]))
{
  printf(STDERR "# total   sequences: 0\n# exit\n");
  exit(0);
}

printf(STDERR "\n# total   sequences: %d taxa: %d\n",
  $#{$cluster_ref}+1,scalar(keys(%intaxa)));
printf(STDERR "# longest sequence: %d (%s)\n",
  $maxlength, (split(/\s+/,$cluster_ref->[$longest_seq][NAME]))[0]);
 

if($#{$cluster_ref} == -1)
{
  warn "# Cannot read input sequences, exit. Please set -P if using a peptide cluster.\n";
  exit;
} 

## read include lists
if($INP_includeA)
{
  # parse include_file A
  open(INCL,$INP_includeA) || die "# EXIT : cannot read $INP_includeA\n";
  while(<INCL>)
  {
    next if(/^#/ || /^$/);
    $taxon = (split)[0];
    $taxon = (split(/\.f/,$taxon))[0];
    $included_input_filesA{$taxon} = 1;
    
    if(!$intaxa{$taxon})
    {
      die "# cannot match $taxon in $INP_clusterfile (included in $INP_includeA)\n";
    }
    else
    {
      foreach $seqid (keys(%id2taxon))
      {
        if($id2taxon{$seqid} eq $taxon)
        {
          push(@listA,$seqid);
        }
      }
    }
  }
  close(INCL); 

  printf(STDERR "# sequences included in group A = %d\n", scalar(@listA));

  # set group 'rest'
  if($INP_includeB)
  {
    # parse include_file B
    open(INCL,$INP_includeB) || die "# EXIT : cannot read $INP_includeB\n";
    while(<INCL>)
    {
      next if(/^#/ || /^$/);
      $taxon = (split)[0];
      $taxon = (split(/\.f/,$taxon))[0];
      $included_input_filesB{$taxon} = 1; 
      
      if(!$intaxa{$taxon})
      {
        die "# cannot match $taxon in $INP_clusterfile (included in $INP_includeB)\n";
      }
      else
      {
        foreach $seqid (keys(%id2taxon))
        {
          if($id2taxon{$seqid} eq $taxon)
          {
            push(@listB,$seqid);
          }
        }
      }
    }
    close(INCL);
  }
  else
  {
    foreach $seqid (keys(%id2taxon))
    {
      next if($included_input_filesA{$id2taxon{$seqid}});
      $included_input_filesB{$taxon} = 1;
      push(@listB,$seqid);
    }
  }
  
  printf(STDERR "# sequences included in group B = %d (rest)\n\n",scalar(@listB));
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
 
if($#{$cluster_ref}+1 < 2)
{
  warn "\n# Need at least two input sequences to be aligned, exit\n";
  exit;
}


## align cluster sequences
# temp files, removed automatically on exit
my ($fhdb,$filenamedb) = tempfile(UNLINK => 1); # blast DB
my ($fhq,$filenameq)   = tempfile(UNLINK => 1); # blast query
my ($fhb,$filenameb)   = tempfile(UNLINK => 1); # blast output
my ($fhr,$filenamer)   = tempfile(UNLINK => 1); # raw alignment file
my ($fha,$filenamea)   = tempfile(UNLINK => 1); # out alignment file
my $mview_seqname;

if($ext_seq ne '')
{
  $fh = $fhq;

  $mview_seqname = $ext_seqname;
  foreach $char (keys(%CHARSFORMVIEW)){
    $mview_seqname =~ s/\Q$char\E/$CHARSFORMVIEW{$char}/g;
  }
  
  print $fh ">$mview_seqname\n$ext_seq\n";    
}

foreach $seq (0 .. $#{$cluster_ref})
{
  if($ext_seq eq '' && $seq == $longest_seq){ $fh = $fhq }
  else{ $fh = $fhdb }
 
  # make sure no forbidden chars are passed to MVIEW
  $mview_seqname = $short_names{$seq};
  foreach $char (keys(%CHARSFORMVIEW)){
    $mview_seqname =~ s/\Q$char\E/$CHARSFORMVIEW{$char}/g;
  }
 
  print $fh ">$mview_seqname\n$cluster_ref->[$seq][SEQ]\n";
}

close($fhb);
close($fhq);
close($fhdb);

if($INP_nucleotides)
{
  executeFORMATDB($filenamedb,1,1); 
  $command = format_BLASTN_command_aligns($filenameq,$filenameb,
    $filenamedb,$DEFEVALUE,$#{$cluster_ref}+1,$DEFBLASTNTASK);
  push(@trash,$filenamedb.'.nsq', $filenamedb.'.nin', $filenamedb.'.nhr'); 
}
else
{
  executeFORMATDB($filenamedb,0,1);
  $command = format_BLASTP_command_aligns($filenameq,$filenameb,$filenamedb,
    $DEFEVALUE,$#{$cluster_ref}+1);
  push(@trash,$filenamedb.'.psq', $filenamedb.'.pin', $filenamedb.'.phr');
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


## format raw alignment  (note that sequence ids are truncated to 60chr, and that MVIEW adds sequence number
$command = "$ENV{'EXE_MVIEW'} -in blast -out fasta $filenameb | perl -lne 's/^>\\d+:/>/; print' > $filenamer"; 
system($command);  

# substitute MVIEW forbidden char place holders and save output in open filehandle $fha
my ($input_seq, @flipped_FASTA);
open(MVIEWRAW,"<",$filenamer) || die "# ERROR: cannot open mview outfile $filenamer\n";
while(my $line = <MVIEWRAW>)
{
  chomp($line);
  if($line =~ /^>/)
  {
    foreach $char (keys(%CHARSFROMMVIEW)){
      $line =~ s/$char/$CHARSFROMMVIEW{$char}/g;
    }

    if($INP_unaligned)
    {
      # find input sequence that matches this header
      $input_seq = -1;
      foreach $seq (0 .. $#{$cluster_ref})
      {
        if($line =~ m/\Q$short_names{$seq}\E/)
	{
	  $input_seq = $seq;
	  last;
	}
      }
      
      if($input_seq == -1)
      {
        warn "# $line not included in output, failed to align to longest/ref sequence\n"
      }
      else 
      {
        # header
	if($do_PFAM)
        {
          push(@flipped_FASTA, ">$short_names{$input_seq} Pfam:$Pfam_full($Pfam_string)")
        }
        else
        {
          push(@flipped_FASTA, ">$short_names{$input_seq}")
        }

        #sequence, flip if required
        $input_seq = $cluster_ref->[$input_seq][SEQ]; 	
        if($line =~ /\+ \-/)
	{ 
	  $input_seq =~ tr/acgtACGT/tgcaTGCA/;
          $input_seq = reverse($input_seq);	
	}
	push(@flipped_FASTA, $input_seq)
      }  
    }

    print $fha "$line\n";

  } else {

    if($INP_unaligned == 0)
    {
      print $fha "$line\n";
    }
  }
}
close(MVIEWRAW);
close($fha);

if($INP_unaligned) {

  # print output unaligned
  if($INP_outfile)
  {
    $fh = FileHandle->new(">$INP_outfile");
    if(!defined($fh))
    {
      die "# cannot create output file $INP_outfile\n";
    }
  }
  else{ $fh = *STDOUT }

  foreach $seq (@flipped_FASTA)
  {
    print $fh "$seq\n"
  }

  if($INP_outfile)
  {
    print STDERR "\n# unaligned (flipped) file: $INP_outfile\n";
    close($fh);
  }

  exit(0);	
}


## extract SNPs, parsimony-informative sites, private variants and trim alignment if required
my ($align_ref,$nSNPS,$npars,$npriv,$SNPs,$pars,$priv,$missA,$missB) = 
  check_variants_FASTA_alignment($filenamea,!$INP_nucleotides,$INP_blunt,\@listA,\@listB);

if(!defined($align_ref->[0]))
{
  printf(STDERR "# aligned sequences: 0\n# exit\n");
  exit(0);
} 
else 
{
  printf(STDERR "# aligned sequences: %d width:   %d\n\n",
    $#{$align_ref}+1,length($align_ref->[0][SEQ])); 
}

if($INP_includeA)
{
  if($INP_nucleotides) 
  {
    printf(STDERR "# alignment sites: SNP=%d parsimony-informative=%d private=%d unaligned A=%d B=%d (%s)\n",
      $nSNPS,$npars,$npriv,$missA,$missB,$INP_clusterfile);
    printf(STDERR "# private sites=%s\n\n",$priv); 
  }
  else
  {
    printf(STDERR "# alignment sites: missense=%d parsimony-informative=%d private=%d unaligned A=%d B=%d (%s)\n",
      $nSNPS,$npars,$npriv,$missA,$missB,$INP_clusterfile);
    printf(STDERR "# private sites=%s\n\n",$priv);
  }
}
else
{
  if($INP_nucleotides)
  {
    printf(STDERR "# alignment sites: SNP=%d parsimony-informative=%d (%s)\n\n",
      $nSNPS,$npars,$INP_clusterfile);
  }
  else 
  {
    printf(STDERR "# alignment sites: missense=%d parsimony-informative=%d (%s)\n\n",
      $nSNPS,$npars,$INP_clusterfile);
  }
}    


## check aligned & unaligned sequences
# 
# This can happen in several cases:
# i) DNA sequences clustered but aminoacid sequences are being annotated (translated fragments shifted)
# ii) longest sequence, which is used to build MSA, does not match some sequences

my ($fhndb,$filenamenotaligndb) = tempfile(UNLINK => 1); # blast DB
my ($fhnq,$filenamenotalignq) = tempfile(UNLINK => 1); # blast output
my ($fhnb,$filenamenotalignb) = tempfile(UNLINK => 1); # out alignment file
my ($alignedOK,$unaligned,%aligned_taxa) = (0,0);

foreach $seq (0 .. $#{$cluster_ref})
{
  $alignedOK = 0;
  foreach $align_seq (0 .. $#{$align_ref})
  {
    if($align_ref->[$align_seq][NAME] =~ m/\Q$short_names{$seq}\E/)
    {
      if($short_names{$seq} =~ m/^.*?\[(\S+?)\]*$/){ $aligned_taxa{$1}++; }
      $alignedOK = 1;
      last;
    }
  }
  
  if($alignedOK == 0)
  {
    $unaligned++;
    print $fhnq ">$short_names{$seq}\n$cluster_ref->[$seq][SEQ]\n";
  }
  else
  {
    print $fhndb ">$short_names{$seq}\n$cluster_ref->[$seq][SEQ]\n";
  }
}    
    
close($fhnq);    
close($fhndb); 

printf(STDERR "# taxa included in alignment: %d\n",scalar(keys(%aligned_taxa)));  
#foreach $taxon (sort(keys(%aligned_taxa))){ warn "## $taxon $aligned_taxa{$taxon}\n"; }

if($unaligned > 0 && $#{$align_ref} > 0)
{
  if($INP_nucleotides)
  {
    executeFORMATDB($filenamenotaligndb,1,1); 
    $command = format_BLASTN_command($filenamenotalignq,$filenamenotalignb,$filenamenotaligndb,
      $DEFEVALUE,5,0,$DEFBLASTNTASK);
      push(@trash,$filenamenotaligndb.'.nsq', $filenamenotaligndb.'.nin', $filenamenotaligndb.'.nhr');
  }
  else
  {
    executeFORMATDB($filenamenotaligndb,0,1); 
    $command = format_BLASTP_command($filenamenotalignq,$filenamenotalignb,$filenamenotaligndb,
      $DEFEVALUE,5,0);
    push(@trash,$filenamenotaligndb.'.psq', $filenamenotaligndb.'.pin', $filenamenotaligndb.'.phr');  
  }

  system($command);
  if($? != 0)
  {
    die "# EXIT: failed while running BLAST search ($command)\n"; 
  }
  else
  {
    warn "\n# sequences not included in multiple alignment (n=$unaligned):\n";
 
    my $unhits = 0;
    my %seen; 
    open(BLASTOUT,'<',$filenamenotalignb) || die "# EXIT: cannot open $filenamenotalignb\n";
    while(<BLASTOUT>)
    {
      my @data = split(/\t/,$_,3);
      warn "# $data[0] (matches $data[1])\n" if(!defined($seen{ $data[0] }));
      $seen{ $data[0] }++;
      $unhits++;
    }
    close(BLASTOUT);

    # show unaligned sequences when there are no double-check hits
    if($unhits == 0)
    {
      open(UNALIGNED,'<',$filenamenotalignq) || die "# EXIT: cannot open $filenamenotalignq\n";
      while(<UNALIGNED>)
      {
        if(/^>(\S+)/){ print "# $1\n"; }
      }
      close(UNALIGNED);  
    }
    
    unlink(@trash);
  }
}

  
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
  print STDERR "\n# alignment file: $INP_outfile\n"; 
  close($fh);
}  


## if requested attempt to collapse overlapping sequences of t
if($INP_collapse && $INP_outfile)
{
  warn "\n# collapsing taxon aligned sequences overlap >= $INP_collapse residues\n";
  my $collapsed_align_ref = 
    collapse_taxon_alignments($INP_outfile,!$INP_nucleotides,$INP_collapse,$MAXMISMCOLLAP,$MAXGAPSCOLLAP);
  
  my $collapsed_outfile_name;
  if($INP_outfile =~ m/(\S+?)\.(\S+)$/)
  {
    $collapsed_outfile_name = $1.'.collapsed.'.$2;
  }
  else{ $collapsed_outfile_name = $INP_outfile.'.collapsed'; }
   
  open(COLLAPSED,">",$collapsed_outfile_name) || die "# cannot create $collapsed_outfile_name\n";
  foreach $seq (0 .. $#{$collapsed_align_ref})
  {
    print COLLAPSED ">$collapsed_align_ref->[$seq][NAME]";
    
    if($do_PFAM){ print COLLAPSED " Pfam:$Pfam_full($Pfam_string)" }
    
    print COLLAPSED "\n$collapsed_align_ref->[$seq][SEQ]\n";
  }
  
  printf( STDERR "\n# collapsed alignment file: %s (aligned sequences: %d)\n",
    $collapsed_outfile_name,$#{$collapsed_align_ref}+1); 
  close(COLLAPSED);
}

