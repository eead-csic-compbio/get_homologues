#!/usr/bin/perl -w
# 2013 Bruno Contreras-Moreira, Pablo Vinuesa

# downloads a list of genomes from the NCBI and its FTP sites
# download_genomes_ncbi.pl

use File::Fetch;
use Net::FTP;
use strict;

my $VERSION  = 1.0;

# http://www.ncbi.nlm.nih.gov/books/NBK25500/
my $EFETCHEXE = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gbwithparts&retmode=text&id=";
my $NCBIHOST  = "ftp.ncbi.nlm.nih.gov";
my $COMPLETE_GENOMES_DIR = "/genomes/genbank/bacteria/";

#my $REFSEQ_GENOMES_DIR   = "/refseq/release/complete/"; # disabled as nucleotide sequences might be omitted
my $WGS_GENOMES_DIR      = "/genbank/wgs/";
my $LOCAL_GENOMES_DIR    ="./";  # where to store downloaded files
my $EXTENSION       = "gbk";
my $REFSEQEXTENSION = "genomic.gbff";
my $WGSEXTENSION    = "gbff";
my $SLEEPTIME = 2;
$File::Fetch::WARN = 0; # no warnings

#######################################################

my ($size,$localsize,$file,$localfile,$listfile,%wanted_genome,%new_name,$newname,$corr_filename);
my ($genome,$dir,$n_of_files,$single_file,$n_of_downloaded_files,$n_of_wanted_files);

if(!$ARGV[0])
{
  die "# usage: $0 v.$VERSION: <genome_list.txt>\n\n".
    "File genome_list.txt must list names of genomes to be downloaded, as explained in the manual. Example:\n".
    "Yersinia_pestis_Pestoides_F Yersinia_pestis_Pestoides # complete genome\n".
    "NC_010159                   Yersinia_pestis_Angola    # known accession number\n".
    "AAYU                        Yersinia_pestis_B42003004 # WGS genome, possibly a draft\n\n"
}
else{ $listfile = $ARGV[0] }

# 0) read list file ###################################
$n_of_wanted_files = 0;
open(LIST,$listfile) || die "# cannot read $listfile\n";
while(<LIST>)
{
  next if(/^#/ || /^\s+/);
  ($genome,$newname) = split;
  $wanted_genome{$genome}{'order'} = $n_of_wanted_files;
  $wanted_genome{$genome}{'wanted'} = 1;
  $wanted_genome{$genome}{'downloaded'} = 0;
  $wanted_genome{$genome}{'appears'}++;
  if($wanted_genome{$genome}{'appears'} > 1)
  {
    die "# $0: error, genome $genome seems to appear twice in input list\n";
  }
  elsif($newname){ $new_name{$genome} = $newname }
  $n_of_wanted_files++;
}
close(LIST);
print "# number of genomes in '$listfile' = $n_of_wanted_files\n";

# 1) look for these genomes ###########################
$n_of_downloaded_files = 0;
chdir($LOCAL_GENOMES_DIR);

# 1.0) first try acession number way
print "\n# checking GenBank accessions with e-utils efetch ...\n";
foreach $genome (sort {$wanted_genome{$a}{'order'}<=>$wanted_genome{$b}{'order'}} keys(%wanted_genome))
{
  $file = $new_name{$genome} || $genome;
  $file = "_$file.$EXTENSION"; #print "$file $wanted_genome{$genome}{'order'}\n";next;
  if(-s $file)
  {
    print ">$file\n";
    $wanted_genome{$genome}{'downloaded'} = 1;
    next;
  }

  eval
  {
    foreach my $g (split(/,/,$genome))
    {
      my ($ff,$gbkcontent);
      $ff = File::Fetch->new(uri =>$EFETCHEXE.$g);
      if($ff->fetch( to => \$gbkcontent ))
      {
        open(GBKFILE,">>$file") || die "# cannot append to $file\n";
        print GBKFILE $gbkcontent;
        close(GBKFILE);
      }
    }

    if(-s $file)
    {
      print ">$file\n";
      if(gbk_contains_CDSs($file))
      {
        $wanted_genome{$genome}{'downloaded'} = 1;
        $n_of_downloaded_files++;
      }
      else
      {
        $wanted_genome{$genome}{'downloaded'} = 1;
        $corr_filename = $file;
        $corr_filename =~ s/^_//;
        rename($file,$corr_filename);
        print "# file $corr_filename does not contain annotated CDSs\n";
      }
    }
    };

  sleep($SLEEPTIME);
}

if($n_of_downloaded_files == $n_of_wanted_files){ die "\n# bye!\n"; }

# 1.1) now browse FTP servers
my $ftp = Net::FTP->new($NCBIHOST,Debug=>0,Timeout=>3600) or die "Cannot connect to $NCBIHOST : $@";

$ftp->login("anonymous",'-anonymous@') or die "# cannot login ", $ftp->message;
$ftp->binary();

# 1.2) see in the complete genomes directory #
print "\n# checking complete genomes...\n";
$ftp->cwd($COMPLETE_GENOMES_DIR) or die "# cannot change working directory ", $ftp->message;

foreach $dir (@{$ftp->dir()})
{
  $dir = (split(/\s+/,$dir))[-1];

  my $diff_name = '';
  foreach $genome (sort {$wanted_genome{$a}{'order'}<=>$wanted_genome{$b}{'order'}} keys(%wanted_genome))
  {
    if($dir =~ /^$genome/)
    {
      $diff_name = $dir;
      last;
    }
  }

  next if(!$wanted_genome{$dir} && $diff_name eq '');

  print ">$dir\n";

  # 1.2.1) cd to this subdirectory
  $ftp->cwd($dir)  or die "# cannot change to subdirectory ", $ftp->message;

  my @files;
  foreach $file (@{$ftp->dir()})
  {
    $file = (split(/\s+/,$file))[-1];
    next if($file !~ /\.$EXTENSION/);  print "$file\n";
    push(@files,$file);

    $size = $ftp->size($file);
    $localfile = $file;

    if(!-s $file)
    {
      print "# copying $file $size ...\n";
      if(!$ftp->get($file)){ die "# Cannot get $file\n" }
    }
    else
    {
      $localsize = (stat($localfile))[7];

      if($size != $localsize)
      {
        print "# updating $file $size ($localsize) ...\n";
        if(!$ftp->get($file)){ die "# Cannot get $file\n" }
      }
      else{ print "# file $file unchanged\n"; }
    }
  }

  # 1.2.2) concatenate files into a single .gbk file
  if(@files)
  {
    if($new_name{$dir})
    {
      $single_file = "_".$new_name{$dir} . ".$EXTENSION";
    }
    else{ $single_file = "_".$dir . ".$EXTENSION"; }

    unlink($single_file) if(-s $single_file);
    foreach $file (@files)
    {
      system("cat $file >> $single_file");
    }

    $ftp->cwd('../') or die "# cannot return to working directory ", $ftp->message;

    if(-s $file && gbk_contains_CDSs($file))
    {
      $wanted_genome{$genome}{'downloaded'} = 1;
      $n_of_downloaded_files++;
    }
    else
    {
      $wanted_genome{$genome}{'downloaded'} = 1;
      $corr_filename = $file;
      $corr_filename =~ s/^_//;
      rename($file,$corr_filename);
      print "# file $corr_filename does not contain annotated CDSs\n";
    }
  }

  sleep($SLEEPTIME);
}

if($n_of_downloaded_files == $n_of_wanted_files){ die "\n# bye!\n"; }

# 1.3) now check the REFSSEQ directory #
# WARNING: disabled as nucleotide sequences might be omitted , Jan2011
#print "\n# checking bacterial REFSEQ genomes...\n";
#$ftp->cwd($REFSEQ_GENOMES_DIR) or die "# Cannot change to REFSEQ directory ", $ftp->message;
#
#foreach $dir (keys(%wanted_genome))
#{
#	next if($wanted_genome{$dir}{'downloaded'});
#
#	my @files;
#	foreach $file (@{$ftp->dir()})
#	{
#		$file = (split(/\s+/,$file))[-1];
#
#		next if($file !~ /$dir/ || $file !~ /\.$REFSEQEXTENSION/); print "$file\n";
#
#		push(@files,$file);
#
#		$size = $ftp->size($file);
#		$localfile = $file;
#
#		if(!-s $file)
#   		{
#      			print "# copying $file $size ...\n";
#      			if(!$ftp->get($file)){ die "# Cannot get $file\n" }
#   		}
#   		else
#   		{
#      			$localsize = (stat($localfile))[7];
#
#      			if($size != $localsize)
#      			{
#         			print "# updating $file $size ($localsize) ...\n";
#         			if(!$ftp->get($file)){ die "# Cannot get $file\n" }
#      			}
#				else{ print "# file $file unchanged\n"; }
#   		}
#	}
#
#	if(@files)
#	{
#		if($new_name{$dir})
#		{
#			$single_file = "_".$new_name{$dir} . ".$EXTENSION";
#		}
#		else{ $single_file = "_".$dir . ".$EXTENSION"; }
#
#		print "> $single_file\n";
#		unlink($single_file) if(-s $single_file);
#		foreach $file (@files)
#		{
#			system("zcat $file >> $single_file");
#		}
#
#		$wanted_genome{$dir}{'downloaded'} = 1;
#		$n_of_downloaded_files++;
#	}
#
#	sleep($SLEEPTIME);
#}
#
#if($n_of_downloaded_files == $n_of_wanted_files){ die "\n# bye!\n"; }

# 1.4) finally check the WGS directory #
print "\n# checking bacterial WGS genomes...\n";
$ftp->cwd($WGS_GENOMES_DIR) or die "# cannot change to WGS directory ", $ftp->message;

foreach $dir (sort {$wanted_genome{$a}{'order'}<=>$wanted_genome{$b}{'order'}} keys(%wanted_genome))
{
  next if($wanted_genome{$dir}{'downloaded'});

  my @files;
  foreach $file (@{$ftp->dir()})
  {
    $file = (split(/\s+/,$file))[-1];

    next if($file !~ /$dir/ || $file !~ /\.$WGSEXTENSION/ || $file =~ /mstr/);
    push(@files,$file); print "$file\n";

    $size = $ftp->size($file);
    $localfile = $file;

    if(!-s $file)
    {
      print "# copying $file $size ...\n";
      if(!$ftp->get($file)){ die "# Cannot get $file\n" }
    }
    else
    {
      $localsize = (stat($localfile))[7];

      if($size != $localsize)
      {
        print "# updating $file $size ($localsize) ...\n";
        if(!$ftp->get($file)){ die "# Cannot get $file\n" }
      }
      else{ print "# file $file unchanged\n"; }
    }
  }

  if(@files)
  {
    if($new_name{$dir})
    {
      $single_file = "_".$new_name{$dir} . ".$EXTENSION";
    }
    else{ $single_file = "_".$dir . ".$EXTENSION"; }

    print "> $single_file\n";
    unlink($single_file) if(-s $single_file);
    foreach $file (@files)
    {
      system("zcat $file >> $single_file");
    }

    if(-s $single_file && gbk_contains_CDSs($single_file))
    {
      $wanted_genome{$dir}{'downloaded'} = 1;
      $n_of_downloaded_files++;
    }
    else
    {
      $wanted_genome{$dir}{'downloaded'} = 1;
      $corr_filename = $single_file;
      $corr_filename =~ s/^_//;
      rename($single_file,$corr_filename);
      print "# file $corr_filename does not contain annotated CDSs\n";
    }
  }

  sleep($SLEEPTIME);
}

print "\n# number of handled genomes = $n_of_downloaded_files (out of $n_of_wanted_files)\n";

$ftp->quit();

sub gbk_contains_CDSs
{
  my ($gbkfile) = @_;
  my $n_of_CDSs = 0;
  open(GBK,$gbkfile) || die "# check_gbk_content: cannot read $gbkfile\n";
  while(<GBK>){ if(/^     CDS /){ $n_of_CDSs++ } }
  close(GBK);
  return $n_of_CDSs;
}
