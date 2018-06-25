#!/usr/bin/env perl

# 2018 Alvaro Rodriguez del Rio and Bruno Contreras-Moreira:
# http://www.eead.csic.es/compbio (Laboratory of Computational and Structural Biology, EEAD-CSIC, Spain)

# This script takes a file with tab-separated BLAST results, with the custom-format used by get_homologues,
# and computes length-normalized E-value & bitscores, as done originally by
# OrthoFinder (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2)

$|=1;

use strict;
use warnings;
use Getopt::Std;
use Benchmark;
use FindBin '$Bin';
use lib $Bin; #"$Bin/../lib";
use Statistics::LineFit; # https://metacpan.org/pod/Statistics::LineFit

my $MINBITSCORE    = 100;  # smaller values are not fitted
my $PERCENTILE     = 95;   # Alvaro's tests suggest 99 might be better for peptides
my $NORMSCALE      = 1000; # normalized bit-scores are multiplied by this constant
my $DBSIZE         = 100_000_000; # used to compute comparable normalized E-values
my $ROUNDEVALUE    = 1e-181; # smaller Evalues set to zero

my (%opts,$blastfile,$outfile,$intercept,$slope,$start,$end);
my ($format,$minbits,$perc,$dbsize,$scale,$round) = 
    ('all',$MINBITSCORE,$PERCENTILE,$DBSIZE,$NORMSCALE,0);

getopts('hri:o:m:p:s:c:f:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-i input BLAST tab-separated file           (required, might be gzipped)\n";
  print   "   expects 12 columns:\n";
  print   "   qseqid sseqid pident length qlen slen\n";
  print   "   qstart qend sstart send evalue bitscore\n";
  print   "-o output filename                          (required)\n";
  print   "-f output format, to choose from:           (optional, default: -f $format)\n";
  print   "   all: all columns\n";
  print   "   norm: evalue & bitscore only\n";
  print   "-m min bitscore for linear model            (optional, default: -m $minbits)\n";
  print   "-p percentile to select top hits in bins    (optional, default: -p $perc)\n";
  print   "-s size of BLAST db in residues             (optional, default: -s $dbsize)\n";
  print   "-c constant to multiply bitscores           (optional, default: -c $scale)\n";
  print   "-r round normalized evalues; if set         (optional, default: no rounding)\n";  
  print   "   evalues<$ROUNDEVALUE are set to 0.0\n";
  exit(0);
}

if(defined($opts{'i'})){  $blastfile = $opts{'i'}; }
else
{ 
  die "# EXIT : need a -i BLAST tab-separated file with 12 cols: \n".
    "# qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore\n"; 
}

if(defined($opts{'o'})){  $outfile = $opts{'o'}; }
else{ die "# EXIT : need output filename (-o)\n"; }

if(defined($opts{'f'}))
{  
  $format = $opts{'f'}; 
  if($format ne 'all' && $format ne 'norm')
  {
    die "# EXIT: accepted formats are: all, norm\n";  
  }
}

if(defined($opts{'m'})){ $minbits = $opts{'m'} }

if(defined($opts{'p'})){ $perc = $opts{'p'} }

if(defined($opts{'s'})){ $dbsize = $opts{'s'} }

if(defined($opts{'c'})){ $scale = $opts{'c'} }

if(defined($opts{'r'})){ $round = $ROUNDEVALUE }

print "# $0 -i $blastfile -f $format -m $minbits -p $perc -s $dbsize -c $scale -r $round\n\n";

$start = new Benchmark();  
 
my ($blast_output_ref, $length_products_ref, $scores_ref ) = 
  parse_pairwise_blast_file($blastfile, $minbits);

my ($top_log10_scores_ref, $top_log10_lengths_ref) = 
  extract_longest_hits($scores_ref, $length_products_ref, $perc);

# fit a linear model of log(bitscore) ~ f(log(length-product))
my $lineFit = Statistics::LineFit->new();
$lineFit->setData($top_log10_lengths_ref,$top_log10_scores_ref) || 
  die "# ERROR: invalid regression data\n";
  
if(defined($lineFit->rSquared()))
{
  ($intercept, $slope) = $lineFit->coefficients();
  printf("# linear model: log(bits) = %1.3f x log(len_prod) + %1.3f\n\n",
    $slope,$intercept);
}
else
{
  die "# ERROR: linear fit failed\n";
}

normalize_scores_to_file($blast_output_ref, $intercept, $slope, 
  $scale, $dbsize, $round, $format, $outfile);

print "# output: $outfile\n";

$end = new Benchmark();  
print "\n# runtime: ".timestr(timediff($end,$start),'all')."\n";
exit(0);



# takes a BLAST tab-separated file with custom-format
# qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
# and parses sequence lengths and raw bit scores. 
# For each pair of sequences Q,S only the first (max) bitscore is taken.
# Returns refs to three arrays:
# 1: reference to array with entries [ bits, Qlength, length-product, raw BLAST string ]
# 2: reference to array with products of sequence length
# 3: reference to array with raw bitscore in same order as 2  
# 
# Note that only hits with bitscore >= $min_bitscore are considered for 2 & 3
#
sub parse_pairwise_blast_file
{
  my ($infile, $min_bitscore) = @_;

  my ($magic,$line,$Qid,$Sid,$Qlength,$Slength,$bits,$lproduct);  # Q=query, S=subject
  my (@blast_output, @length_products, @scores, %BLASTDB);
  
  # check input file format and open it accordingly
  open(INFILE,$infile) || die "# parse_pairwise_blast_file: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    if(!open(INFILE,"gzip -dc $infile |"))
    {
      die "# parse_pairwise_blast_file: cannot read GZIP compressed $infile $!\n"
        ."# please check gzip is installed\n";
    }
  }
  else{ open(INFILE,"<$infile") || die "# parse_pairwise_blast_file: cannot read $infile $!\n"; }
    
  while($line = <INFILE>)
  {
    chomp($line);
    
    #qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
    if($line =~ m/(\S+)\s(\S+)\s\S+\s\S+\s(\S+)\s(\S+)\s\S+\s\S+\s\S+\s\S+\s\S+\s(\S+)/){

      ($Qid,$Sid,$Qlength,$Slength,$bits) = ($1,$2,$3,$4,$5);

      $lproduct = $Qlength * $Slength;

      # store relevant elems & original BLAST output
      push(@blast_output, [ $bits, $Qlength, $lproduct, $line ]);

      next if($bits < $min_bitscore || defined($BLASTDB{$Qid}{$Sid})); # skip poor/secondary matches
      
      push (@length_products, $lproduct);
      push (@scores, $bits);
      
      # remember this Q,S pair
      $BLASTDB{$Qid}{$Sid} = 1;
    }
  }
  
  close(INFILE);
  
  return (\@blast_output, \@length_products, \@scores);
}

# Extracts length products and scores of top/longest BLAST hits. 
# Contains hard-coded cut-offs for making bins of BLAST hits.
# Takes three params:
# 1: ref to array with bitscores
# 2: ref to array with sequence length products matching bitscores in 1
# 3: percentile to extract
# Returns:
# 1: ref to array with log10(bitscores) of top hits
# 2: ref to array with log10(length products) of top hits
sub extract_longest_hits
{
  my ($scores_ref, $length_ref, $topPer) = @_;
  
  my ($hit,$bits,$length,$nInBins,@topScores,@topLengths);
  my $nBLASThits = scalar(@$scores_ref);
  my $last_hit = $nBLASThits-1;

  sub log10 
  {
    my $n = $_[0];
    return log($n)/log(10);
  }
  
  # get sorted index (by length product) for later use
  my @sortedIndex = sort { $length_ref->[$a] <=> $length_ref->[$b] } 0 .. $last_hit;

  # calculate number of blast hits per bin
  if($nBLASThits<100)
  {
    # if too few take them all as top hits and exit
    foreach $hit (@sortedIndex)
    {
      push(@topScores, log10($scores_ref->[$hit]));
      push(@topLengths, log10($length_ref->[$hit]));
    }  
  
    return (\@topScores, \@topLengths);
  }
  elsif($nBLASThits > 5000)
  {
    $nInBins = 1000;
  }
  elsif($nBLASThits > 1000)
  {
    $nInBins = 200;
  }
  else
  {
    $nInBins = 20;
  }
  
  # extract top hits in each bin
  my ($first,$last);
  my $nBins = int($nBLASThits/$nInBins)+1;
  my $hitsToExtract = int(((100-$topPer)/100)*$nInBins);
  
  for(my $i = 0; $i < $nBins; $i++)
  {
    my (@theseLengths,@theseScores,@cutoffScores,@cutoffLengths);
  
    # get coordinates of the elements which will go into the bin
    $first = $i*$nInBins;
    $last = $first + $nInBins - 1; 

    # check size of bin
    if($i == $nBins-1) # last bin
    { 
      $last = $last_hit;

      if($last <= $first) 
      {
        # no more bins, exit
        return (\@topScores, \@topLengths);
      }
      
      # recalculate the number of values to extract from the last bin, 
      # which will be lower than the normal size; min=1
      $hitsToExtract = int(((100-$topPer)/100) * ($last_hit-$first)) || 1;
    } #print "$i $first $last $last_hit $hitsToExtract\n";

    # sort data of this bin by bit score
    my @sortedIndexBin = sort { $scores_ref->[$b] <=> $scores_ref->[$a] } @sortedIndex[ $first .. $last ];

    # extract top scores of this bin
    foreach $hit (@sortedIndexBin[ 0 .. $hitsToExtract-1])
    {
      push(@topScores, log10($scores_ref->[$hit]));
      push(@topLengths, log10($length_ref->[$hit]));
    }
  }
  
  return (\@topScores, \@topLengths);
}


# Loops through the list of BLAST hits in original order and 
# prints the normalized scores in the requested format 
# Takes 9 params:
# 1: ref to array with entries [ bits, Qlength, length-product, raw BLAST string ]
# 2: intercept of linear model
# 3: slope of linear model relating bitscore to length-product
# 4: scale of normalized bitscore
# 5: size of BLAST db in residues
# 6: cutoff to round E-values
# 7: output format
# 8: filename for output
sub normalize_scores_to_file
{
  my ( $blast_output_ref, $intercept, $slope, $scale, $dbsize, $round_cutoff, $format, $outfilename ) = @_;

  my ($norm_bits,$norm_evalue);

  open(OUTFILE,">",$outfilename) || 
    die "# normalize_scoreis_to_file: cannot create $outfilename\n";

  foreach my $row (@$blast_output_ref)
  {
    # $row is a ref to a vector of 4 elems: bits, Qlength, lproduct, BLAST_tab_string
    
    $norm_bits = sprintf("%4.0f", $scale * $row->[0] / ((10**$intercept)*( $row->[2] ** $slope)));

    # compute E-value from normalized bits (https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html)
    $norm_evalue = sprintf("%1.3g",($row->[1]*$dbsize)/(2**$norm_bits));

    if($norm_evalue < $round_cutoff){ $norm_evalue = 0.0 }
    
    if($format eq 'norm')
    {
      print OUTFILE "$norm_evalue\t$norm_bits\n";
    }
    else
    {
      $row->[3] =~ s/\t\S+\t\S+$//; # delete original E-value & bitscore
      print OUTFILE "$row->[3]\t$norm_evalue\t$norm_bits\n";
    }
  }
  
  close(OUTFILE);
}
