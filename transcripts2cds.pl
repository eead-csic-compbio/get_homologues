#!/usr/bin/env perl

# 2015-6 Bruno Contreras-Moreira (1)
# 1: http://www.eead.csic.es/compbio (Laboratory of Computational Biology, EEAD/CSIC, Spain)

# This script can be used to extract coding sequences encoded in input transcripts

# Uses some code from third-party software:
# TransDecoder (http://transdecoder.sourceforge.net), primarily maintained by Brian Haas at the Broad Institute and
# Alexie Papanicolaou at the Commonwealth Scientific and Industrial Research Organisation (CSIRO).

$|=1;

use strict;
use Getopt::Std;
use File::Basename;
use Cwd;
use FindBin '$Bin';
use lib "$Bin/lib";
use lib "$Bin/lib/est";
use lib "$Bin/lib/est/TransDecoder_r20140704/PerlLib";
use lib "$Bin/lib/bioperl-1.5.2_102";
use phyTools;
use transcripts;
#use cppcode; #CPP
use Bio::Tools::CodonTable;
use Nuc_translator;

my $VERSION = 1.1;

my @FEATURES2CHECK = ('EXE_BLASTX_EST','EXE_FORMATDB_EST','EXE_TRANSDECOD_EST'); #,'EXE_GMAP','EXE_GMAPBUILD');

my ($n_of_cpus,$INP_minORFlength) = ($BLASTX_NOCPU,$MINORFLENGTH);
my ($evalue_cutoff,$INP_blastdb) = ($MAXBLASTXEVALUE,$BLASTXDB);
my ($INP_plus_strand,$INP_gencode,$input_FASTA_file,$input_reference_FASTA_file) = (0,1);
my ($seq_transcod,$seq_transcod_prot,$seq_blastx,$seq_cds,$seq_prot,$besthit);
my ($ref_blastxseqs_cds,$ref_blastxhits,$ref_transcodseqs_cds,$ref_transcodseqs_prot);
my ($ref_gmap_cds,$ref_gmap_besthit,$gmap_besthit,$intronless_infile);
my ($root,$tmproot,$command,$output_mask,$output_mask2,$output_mask3); #mask2=transdecoder,3=blastx
my ($n_of_ORFs,$n_of_noORFs,$evidence,$seq,$seqname,@input_files,%opts);

getopts('hvGpn:l:d:E:g:m:', \%opts);

if(($opts{'h'})||(!$opts{'G'} && scalar(@ARGV)==0))
{
  print   "\nusage: $0 [options] <input FASTA file(s) with transcript nucleotide sequences>\n\n";
  print   "-h this message\n";
  print   "-p check only 'plus' strand                                  (optional, default both strands)\n";
  print   "-l min length for CDS                                        (optional, default=$INP_minORFlength amino acid residues)\n";
  print   "-g genetic code to use during translation                    (optional, default=$INP_gencode, example: -g 11)\n";
  #print   "-m map transcripts against reference FASTA genomic file      (optional, gets common genomic coordinates and splices introns)\n";
  print   "-d run blastx against selected protein FASTA database file   (default=swissprot, example: -d /path/to/sequences.faa)\n";
  print   "-E max E-value during blastx search                          (default=$evalue_cutoff)\n";
  print   "-n number of threads for BLASTX jobs                         (default=$n_of_cpus)\n\n";
  print   "-G show available genetic codes and exit\n";
  exit(0);
}

if(defined($opts{'v'}))
{
  print "\n$0 version $VERSION (2016)\n";
  print "\nProgram written by Bruno Contreras-Moreira\n";
  print "\nhttp://www.eead.csic.es/compbio (Estacion Experimental Aula Dei/CSIC/Fundacion ARAID, Spain)\n";
  print "\nThis software employs binaries from different authors, please cite them accordingly:\n";
  print " NCBI Blast-2.2 (blast.ncbi.nlm.nih.gov , PubMed=9254694,20003500)\n";
  print " TransDecoder (http://transdecoder.sourceforge.net)\n";
  print " Bioperl v1.5.2 (www.bioperl.org , PubMed=12368254)\n";

  # check all binaries and data needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'G'}))
{
  my $codontables = Bio::Tools::CodonTable->tables();
  foreach my $id (sort {$a<=>$b} keys(%$codontables))
  {
    print "$id\t:\t$codontables->{$id}\n";
  }
  exit(0);
}

# check selected genetic code
if(defined($opts{'g'}))
{
  $INP_gencode = $opts{'g'};

  my $codontables = Bio::Tools::CodonTable->tables();
  if(!$codontables->{$INP_gencode})
  {
    print "# EXIT : need a valid genetic code table, please choose choose one from the list:\n";
    foreach my $id (sort {$a<=>$b} keys(%$codontables))
    {
      print "$id\t:\t$codontables->{$id}\n";
    }
    exit(-1);
  }
  else
  {
    # set genetic code table for translations
    Nuc_translator::use_specified_genetic_code($INP_gencode);
    $output_mask  .= "_gencode$INP_gencode";
    $output_mask2 .= "_gencode$INP_gencode";
    $output_mask3 .= "_gencode$INP_gencode";
  }
}

# check input files
foreach $input_FASTA_file (@ARGV)
{
  if($input_FASTA_file =~ /([\+])/)
  {
    die "# EXIT : need a valid input FASTA file, offending char: '$1'\n";
  }
  elsif(!-e $input_FASTA_file)
  {
    die "# EXIT : need a valid input FASTA file\n";
  }
  else{ push(@input_files,$input_FASTA_file) }
}

# check rest of parameters
if(defined($opts{'n'}) && $opts{'n'} > 0)
{
  $n_of_cpus = $opts{'n'};
  $BLASTX_NOCPU = $n_of_cpus;
}

if(defined($opts{'p'}))
{
  $INP_plus_strand = 1;
  $output_mask  .= '_str1';
  $output_mask2 .= '_str1';
  $output_mask3 .= '_str2';
}

if(defined($opts{'l'}) && $opts{'l'} > 0)
{
  $INP_minORFlength = $opts{'l'};
  $MINORFLENGTH = $INP_minORFlength;
}
$output_mask .= "_l$INP_minORFlength";
$output_mask2 .= "_l$INP_minORFlength";

if(defined($opts{'E'}) && $opts{'E'} > 0)
{
  $evalue_cutoff = $opts{'E'};
  $MAXBLASTXEVALUE = $evalue_cutoff;
}
$output_mask  .= "_E$evalue_cutoff";
$output_mask3 .= "_E$evalue_cutoff";

if(defined($opts{'d'}) && glob("$opts{'d'}*"))
{
  $INP_blastdb = $opts{'d'};
  $output_mask  .= '_'.basename($INP_blastdb);
  $output_mask3 .= '_'.basename($INP_blastdb);
}

if(defined($opts{'m'}) && -s $opts{'m'})
{
  $input_reference_FASTA_file = $opts{'m'};
  $output_mask2 .= "_ref$input_reference_FASTA_file";
  $output_mask  .= "_ref$input_reference_FASTA_file";
}

print "# $0 -p $INP_plus_strand -m $input_reference_FASTA_file -d $INP_blastdb ".
  "-E $evalue_cutoff -l $INP_minORFlength -g $INP_gencode -n $n_of_cpus\n";

if(@input_files)
{
  print "# input files(s):\n";
  foreach $input_FASTA_file (@input_files){ print "# $input_FASTA_file\n" }
}

##########################################################################

foreach $input_FASTA_file (@input_files)
{
  my ($compressOK,$tmp_FASTA_file,@trash) = (0);

  print "\n## processing file $input_FASTA_file ...\n";

  $root = basename($input_FASTA_file);
  $tmproot = '_'.$root;

  ## 0) make sure headers of input file are suitable
  my $short_header_file = $tmproot.'.short';
  if(!shorten_headers_FASTA_file($input_FASTA_file,$short_header_file))
  {
    print "# cannot read $input_FASTA_file, skip it\n";
    next;
  }
  push(@trash,$short_header_file);
  
  # from now on use short temp file
  $input_FASTA_file = $short_header_file;
  $compressOK = 0;

  ## 1) run transdecoder on input sequences
  my $transdecod_outfile_prot = $tmproot.$output_mask2.'.transdecoder.pep.gz';
  my $transdecod_outfile_cds  = $tmproot.$output_mask2.'.transdecoder.cds.gz';

  ## optionally run gmap to extract CDSs from transcripts mapped on a genomic reference
  if($input_reference_FASTA_file)
  {
    my $gmaproot = '_'.basename($input_reference_FASTA_file);
    my $tmp_gmapfile = $tmproot.$gmaproot.'.gmap';
    my $gmap_outfile = $tmp_gmapfile.'.gz';
    my $gmap_cds_outfile = $tmproot.$gmaproot.'.gmap.cds.fna';
    $intronless_infile = $gmap_cds_outfile.'.gz';

    if(!-e $intronless_infile)
    {
      if(!-s $gmaproot.".gmapdb/$gmaproot.gmapdb.sarray")
      {
        executeGMAPBUILD('./',"$gmaproot.gmapdb",$input_reference_FASTA_file);
        if(!-s $gmaproot.".gmapdb/$gmaproot.gmapdb.sarray")
        {
          die "# EXIT: cannot format GMAP sequence base $input_reference_FASTA_file\n"
        }
      }

      if(!-e $gmap_outfile)
      {
        if($compressOK)
        {
          $tmp_FASTA_file = extract_compressed_file($input_FASTA_file,$compressOK);
          push(@trash,$tmp_FASTA_file);
        }
        else{ $tmp_FASTA_file = $input_FASTA_file }

        $command = format_GMAP_command('./',$gmaproot.'.gmapdb',$tmp_FASTA_file,$tmp_gmapfile,$n_of_cpus,1);
        print("\n# running gmap...\n");
        system("$command");
        if($? != 0)
        {
          die "# EXIT: failed while running GMAP ($command)\n";
        }

        print "# parsing GMAP output ($tmp_gmapfile) ...\n";
        ($ref_gmap_cds,$ref_gmap_besthit) = parse_GMAP_results($tmp_FASTA_file,$tmp_gmapfile,1);

        system("gzip $tmp_gmapfile");
      }
      else
      {
        print "# parsing re-used GMAP output ($gmap_outfile) ...\n";
        ($ref_gmap_cds,$ref_gmap_besthit) = parse_GMAP_results($input_FASTA_file,$gmap_outfile,1);
      }

      # create file with intronless sequences for transdecoder
      open(INTRONLESS,">$gmap_cds_outfile") || die "# EXIT: cannot create $gmap_cds_outfile\n";
      foreach $seqname (keys(%$ref_gmap_cds))
      {
        if($ref_gmap_besthit->{$seqname}){ $gmap_besthit = 'genomic_match:'.$ref_gmap_besthit->{$seqname} } else { $gmap_besthit = '' }
        printf(INTRONLESS ">%s %s\n%s\n",$seqname,$gmap_besthit,$ref_gmap_cds->{$seqname});
      }
      close(INTRONLESS);

      system("gzip $gmap_cds_outfile");
    }
    else
    {
      print "# re-using intronless sequences in $intronless_infile ...\n";

      my $fasta_ref = read_FASTA_file_array($intronless_infile);
      foreach $seq ( 0 .. $#{$fasta_ref} )
      {
        if($fasta_ref->[$seq][NAME] =~ m/genomic_match:(\S+)/)
        {
          $seqname = (split(/\s+/,$fasta_ref->[$seq][NAME]))[0];
          $ref_gmap_besthit->{$seqname} = $1;
        }
      }
    }
  }

  # now perform a plain transdecoder run
  if(!-e $transdecod_outfile_cds)
  {
    my $tmp_FASTA_file_trans;

    if($intronless_infile) # files with ref intron-filtered sequences, gzipped by design
    {
      $tmp_FASTA_file_trans = extract_compressed_file($intronless_infile,'gzip');
      push(@trash,$tmp_FASTA_file_trans);
    }
    else
    {
      if($compressOK)
      {
        $tmp_FASTA_file = extract_compressed_file($input_FASTA_file,$compressOK);
        push(@trash,$tmp_FASTA_file);
      }
      else{ $tmp_FASTA_file = $input_FASTA_file }

      $tmp_FASTA_file_trans = $tmp_FASTA_file;
    }

    print "# running transdecoder...\n";

    $command = format_transdecoder_command($tmp_FASTA_file_trans,$INP_minORFlength,$INP_plus_strand,$INP_gencode);
    system($command);
    if($? != 0)
    {
      die "# ERROR: failed while running transdecoder ($command)\n";
    }
    else
    {
      # compress and rename outfiles with mask2
      my $tmpcdsfile = basename($tmp_FASTA_file_trans).'.transdecoder.cds';
      my $tmppepfile = basename($tmp_FASTA_file_trans).'.transdecoder.pep';
      system("gzip $tmpcdsfile $tmppepfile");
      rename($tmppepfile.'.gz',$transdecod_outfile_prot);
      rename($tmpcdsfile.'.gz',$transdecod_outfile_cds);

      print "# parsing transdecoder output ($transdecod_outfile_cds) ...\n";
      $ref_transcodseqs_cds = parse_transcoder_sequences($transdecod_outfile_cds);
      $ref_transcodseqs_prot = parse_transcoder_sequences($transdecod_outfile_prot);
    }
  }
  else
  {
    print "# parsing re-used transdecoder output ($transdecod_outfile_cds) ...\n";
    $ref_transcodseqs_cds = parse_transcoder_sequences($transdecod_outfile_cds);
    $ref_transcodseqs_prot = parse_transcoder_sequences($transdecod_outfile_prot);
  }

  ## 2) now run blastx to compare transcripts to a large protein set
  if($INP_blastdb)
  {
    my $blastx_outfile    = $tmproot.$output_mask3.'.blastx';
    my $blastx_outfile_gz = $blastx_outfile.'.gz';

    # format sequence database if required
    if(!-s $INP_blastdb .".pal")
    {
      executeFORMATDB_EST($INP_blastdb);
    }

    if(!-s $INP_blastdb .".pal")
    {
      die "# WARNING: cannot format sequence database $INP_blastdb\n";
    }
    else
    {
      if(!$tmp_FASTA_file)
      {
        if($compressOK)
        {
          $tmp_FASTA_file = extract_compressed_file($input_FASTA_file,$compressOK);
          push(@trash,$tmp_FASTA_file);
        }
        else{ $tmp_FASTA_file = $input_FASTA_file }
      }

      if(!-s $blastx_outfile_gz)
      {
        print "# running blastx...\n";
        $command = format_BLASTX_command($tmp_FASTA_file,$blastx_outfile,$INP_blastdb,$evalue_cutoff,$INP_gencode,$INP_plus_strand);
        system($command);
        if($? != 0)
        {
          die "# ERROR: failed while running BLASTX ($command)\n";
        }
        else
        {
          system("gzip $blastx_outfile");
          print "# parsing blastx output ($blastx_outfile_gz) ...\n";
          ($ref_blastxseqs_cds,$ref_blastxhits) =
            parse_blastx_cds_sequences($tmp_FASTA_file,$blastx_outfile_gz);
        }
      }
      else
      {
        print "# parsing re-used blastx output ($blastx_outfile_gz) ...\n";
        ($ref_blastxseqs_cds,$ref_blastxhits) =
          parse_blastx_cds_sequences($tmp_FASTA_file,$blastx_outfile_gz);
      }
    }
  }

  ## 3) output (partial) sequences for each ORF and also transcripts with no ORFs
  my $outfile_cds   = $root.$output_mask.'.cds.fna';
  my $outfile_prot  = $root.$output_mask.'.cds.faa';
  my $outfile_rna   = $root.$output_mask.'.transcript.fna';
  my $outfile_noORF = $root.$output_mask.'.noORF.fna';

  open(FAA,">$outfile_prot") || die "# ERROR: cannot write to $outfile_prot\n";
  open(CDS,">$outfile_cds")  || die "# ERROR: cannot write to $outfile_cds\n";
  open(RNA,">$outfile_rna")  || die "# ERROR: cannot write to $outfile_rna\n";
  open(NORF,">$outfile_noORF")  || die "# ERROR: cannot write to $outfile_noORF\n";

  print "# calculating consensus sequences ...\n";

  #my $consensus = new cppcode::Consensus();         #CPP
  #$consensus->set_priority(2);                      #CPP
  #$consensus->set_min_overlap($MINCONOVERLAP);      #CPP
  #$consensus->set_sources('transdecoder','blastx'); #CPP

  my $fasta_ref = read_FASTA_file_array($tmp_FASTA_file);
  $n_of_ORFs = $n_of_noORFs = 0;
  foreach $seq ( 0 .. $#{$fasta_ref} )
  {

    #$seqname = $fasta_ref->[$seq][NAME]; #print ">$seqname\n";
    $seqname = (split(/\s+/,$fasta_ref->[$seq][NAME]))[0];

    $seq_transcod      = $ref_transcodseqs_cds->{$seqname} || '';
    $seq_transcod_prot = $ref_transcodseqs_prot->{$seqname} || '';
    $seq_blastx        = $ref_blastxseqs_cds->{$seqname} || '';
    $seq_cds = '';

    if($ref_blastxhits->{$seqname}){ $besthit = 'match:'.$ref_blastxhits->{$seqname} } else { $besthit = '' }
    if($ref_gmap_besthit->{$seqname}){ $gmap_besthit = 'genomic_match:'.$ref_gmap_besthit->{$seqname} } else { $gmap_besthit = '' }

    # reconciliate transdecod $ blastx cds sequences if possible
    if($seq_transcod ne '' && $seq_blastx ne '')
    {
      if(substr($seq_transcod_prot,-1) eq '*')
      {
        $seq_transcod = substr($seq_transcod,0,-3);
      }

      ($seq_cds,$evidence) = sequence_consensus($seq_transcod,$seq_blastx,'transdecoder','blastx',2); # Perl

      #$consensus->set_sequences($seqname,$seq_transcod,$seq_blastx); #CPP
      #$consensus->calc_consensus();                                  #CPP
      #$seq_cds  = $consensus->get_consensus();                       #CPP
      #$evidence = $consensus->get_evidence();                        #CPP
    }
    elsif($seq_transcod ne '')
    {
      $seq_cds = $seq_transcod;
      $evidence = 'transdecoder';
    }
    elsif($seq_blastx ne '')
    {
      $seq_cds = $seq_blastx;
      $evidence = 'blastx';
    }

    if((length($seq_cds)/3) < $INP_minORFlength)
    {
      $n_of_noORFs++;
      print NORF ">$seqname\n".$fasta_ref->[$seq][SEQ]."\n";
      next;
    }

    $seq_prot = translate_sequence($seq_cds,1);

    print FAA ">$seqname evidence:$evidence $besthit $gmap_besthit\n$seq_prot\n";
    print CDS ">$seqname evidence:$evidence $besthit $gmap_besthit\n$seq_cds\n";
    print RNA ">$seqname evidence:$evidence $besthit $gmap_besthit\n".$fasta_ref->[$seq][SEQ]."\n";

    $n_of_ORFs++;
  }

  close(FAA);
  close(CDS);
  close(RNA);
  close(NORF);

  printf("# input transcripts = %d\n",scalar(@$fasta_ref));
  printf("# transcripts with ORFs = %d\n",$n_of_ORFs);
  printf("# transcripts with no ORFs = %d\n",$n_of_noORFs);
  print "# output files: ";

  if($n_of_ORFs)
  {
    printf("%s , %s , %s ",$outfile_rna,$outfile_cds,$outfile_prot);
  }
  if($n_of_noORFs)
  {
    printf(", %s",$outfile_noORF);
  }
  else{ unlink($outfile_noORF) }

  print "\n";

  ## 5) clean any tmp files
  unlink(@trash);

} ## foreach
