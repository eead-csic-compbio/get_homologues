#!/usr/bin/env perl

# Script that checks/compiles software required by get_homologues[-est] and 
# checks dependencies for first-time users.
# last checked Aug2024

use strict;
use warnings;
use Cwd;
use FindBin '$Bin';
use File::Path qw(remove_tree);
use lib "$Bin/lib";
use lib "$Bin/lib/est";
use lib "$Bin/lib/bioperl-1.5.2_102/";
require phyTools;
use transcripts;

$|=1;

my $DOWNLOADEXE = 'wget'; # add path if required, curl in MacOS
my $BINTGZFILE = 'bin.tgz';
my $BINOSXTGZFILE = 'binosx.tgz';
my $BINURL = "https://github.com/eead-csic-compbio/get_homologues/releases/download/v3.7/$BINTGZFILE";
my $BINOSXURL = "https://github.com/eead-csic-compbio/get_homologues/releases/download/v3.7/$BINOSXTGZFILE";

my $PFAMSERVERURL   = 'ftp.ebi.ac.uk';
my $PFAMFOLDER      = 'pub/databases/Pfam/current_release/';
my $PFAMHMMFILE     = 'Pfam-A.hmm.gz';
my $BLASTSERVERPATH = 'ftp.ncbi.nih.gov';
my $BLASTEXEPATH    = $BLASTSERVERPATH.'/blast/executables/blast+/LATEST/';
my $SWISSPROTSERVER = 'ftp.uniprot.org';
my $SWISSPROTFOLDER = '/pub/databases/uniprot/current_release/knowledgebase/complete/';
my $SWISSPROTFILE   = $ENV{'BLASTXDB'}.'.gz';

my @supportedOS = ('debian','redhat','suse','macOSX');
my %manager = ('debian'=>'apt-get -y ','redhat'=>'yum -y ','suse'=>'zypper --assume-yes ','macOSX'=>'fink ');
my %packages =
  (
  'R'=>['r-base r-base-dev','R','R-patched R-patched-devel','r-base'],
  'GD'=>['libgd-gd2-perl','perl-GD','perl-GD','gd']
  );

my $nonCOGS = 'diamond-0.8.25 hmmer-3.1b2 mcl-14-137 ncbi-blast-2.14.0+ phylip-3.695';

my ($SOguess,$output,$command,$cwd) = ('','','',getcwd());
my ($downloadOK,$force_unsupervised,$noDBs) = (0,0,0);
my ($COGSonly) = (0);

# unsupervised, forced install
if(defined($ARGV[0]))
{ 
  if($ARGV[0] eq 'force')
  {
    $force_unsupervised = 1;   
    print "# \$force_unsupervised=$force_unsupervised\n\n";
  }
  elsif($ARGV[0] eq 'no_databases')
  { 
    $noDBs = 2;   
    print "# \$no_databases=$noDBs\n\n";
  }
  elsif($ARGV[0] eq 'swissprot')
  {
    $noDBs = 1;
    print "# \$no_databases=$noDBs (swissprot)\n\n";
  }
  elsif($ARGV[0] eq 'test')
  {
    print "# testing only\n";
    exit(0);
  }
} 

if(defined($ARGV[1]))
{
  if($ARGV[1] eq 'COGS')
  {
    $COGSonly = 1
  }
}

##############################################################################################

print "\n### 1) checking required parts: \n\n";

# find out OS
if($^O =~ /darwin/){ $SOguess = $#supportedOS }
else
{
  foreach my $os (0 .. $#supportedOS)
  {
    $output = `$manager{$supportedOS[$os]} 2>&1`; 
    if($output && $output !~ /command not found/)
    {
      $SOguess = $os;  
      last;
    }
  }
}

# check whether bin dir is empty
print "## checking whether source and binaries of dependencies are in place\n";
if(! -e $ENV{'MARFIL'}.'/bin/COGsoft/' || !-e $ENV{'EXE_BLASTP'})
{
  my $userword = '';
  if($noDBs==0 && $force_unsupervised==0)
  {
    print "# cannot locate source and binaries, would you like to download it now? [Y/n]\n";
    $userword = <STDIN>;
  }    

  if($noDBs || $force_unsupervised || $userword =~ m/Y/i )
  {
    chdir($ENV{'MARFIL'}.'/bin/');
    my ($ff,$fpath);
  
    # check OS
    if($supportedOS[$SOguess] eq 'macOSX')
    {
       $BINTGZFILE =  $BINOSXTGZFILE;
       $BINURL = $BINOSXURL;

       # check curl is available
       if(`which curl`)
       {
         $DOWNLOADEXE = 'curl -LJO '
       } 
       elsif(!`which $DOWNLOADEXE`) # check wget is available instead
       {
         die "# need wget, install it with 'brew install wget' and re-run\n";
       }
    }

    if(! -s $BINTGZFILE)
    {
      print "# retrieving $BINURL ...\n";
      system("$DOWNLOADEXE $BINURL");
      if(! -s $BINTGZFILE)
      {
        die "# cannot download file $BINTGZFILE \n\n".
              "<< You might manually download it from $BINURL to $ENV{'MARFIL'}/bin/\n".
              "<< and re-run\n";
      } 
    }
   
    if(-s $BINTGZFILE)    
    {
      print "# extracting $BINTGZFILE ...\n";
      system("tar xfz $BINTGZFILE");

      if($COGSonly)
      {
        mkdir('remove/');
        system("mv $nonCOGS remove/");	
        remove_tree('remove/');
      }

      unlink($BINTGZFILE);
      chdir($cwd);
      $downloadOK = 1;
    }         
  }          
}
else{ print ">> OK\n"; }


# check required precompiled binaries

if(!$COGSonly)
{
  print "\n## checking mcl (lib/phyTools: \$ENV{'EXE_MCL'})\n";

  if($downloadOK){
    $ENV{'EXE_MCL'} = $ENV{'MARFIL'}.'/bin/mcl-14-137/src/shmcl/mcl'
  }

  $output = `$ENV{'EXE_MCL'} 2>&1`;
  if(!$output || $output !~ /usage/)
  {
    print "\n# compiling mcl-14-137 ...\n"; # requires gcc

    if(!-s $ENV{'MARFIL'}.'/bin/mcl-14-137')
    {
      warn "# mcl binaries are missing,\n".
        "# please re-run perl install.pl and type Y when asked to download source and binaries\n";
      exit(-1);   
    }
  
    chdir($ENV{'MARFIL'}.'/bin/mcl-14-137');
    system("./configure 2>&1");
    system("make clean 2>&1");
    system("make 2>&1");
    chdir($cwd); 

    $output = `$ENV{'EXE_MCL'} 2>&1`;
    if($output !~ /usage/)
    {
      die "<< Cannot run $ENV{'EXE_MCL'}, please check the manual for compilation instructions and re-run.\n";
    }
    else{ print ">> OK\n"; }
  }
  else{ print ">> OK\n"; }
}

my %cogs = ("EXE_COGLSE"=>'COGlse',"EXE_MAKEHASH"=>'COGmakehash',
  "EXE_READBLAST"=>'COGreadblast',"EXE_COGTRI"=>'COGtriangles');

foreach my $exe (keys(%cogs))
{
  print "## checking COGsoft/$cogs{$exe} (lib/phyTools: \$ENV{'$exe'})\n";

  if($downloadOK){
    $ENV{$exe} = $ENV{'MARFIL'}."/bin/COGsoft/$cogs{$exe}/$cogs{$exe}"
  }

  $output = `$ENV{$exe} 2>&1`;
  if(!$output || ($output !~ /Usage/ && $output !~ /Error/ && $output !~ /open output/))
  {
    print "\n# compiling COGsoft/$cogs{$exe} ...\n"; # requires g++

    if(!-s $ENV{'MARFIL'}.'/bin/COGsoft/'.$cogs{$exe})
    {
      warn "# COGS binaries are missing,\n".
      "# please re-run perl install.pl and type Y when asked to download source and binaries\n";
      exit(-1);
    }

    chdir($ENV{'MARFIL'}.'/bin/COGsoft/'.$cogs{$exe});
    system("make clean; make 2>&1");
    chdir($cwd);

    $output = `$ENV{$exe} 2>&1`;
    if(!$output || ($output !~ /Usage/ && $output !~ /Error/ && $output !~ /open output/))
    {
      die "<< Cannot run $ENV{$exe}, please check the manual for compilation instructions and re-run ($output).\n";
    }
  }
  else{ print ">> OK\n"; }
}

if(!$COGSonly)
{
  # NCBI blast+
  print "## Checking blast (lib/phyTools: \$ENV{'EXE_BLASTP'})\n";

  if($downloadOK){
      $ENV{'EXE_BLASTP'} = $ENV{'MARFIL'}."/bin/ncbi-blast-2.14.0+/bin/blastp"
  }

  $output = `$ENV{'EXE_BLASTP'} 2>&1`;
  if(!$output || $output !~ /BLAST/)
  {
    die "<< Cannot run shipped blastp, please download it from $BLASTEXEPATH ,\n".
      "<< install it and edit variable BLAST_PATH as explained in the manual\n".
      "<< (inside set_phyTools_env in file lib/phyTools.pm) .\n".
      "<< Then re-run\n";
  }
  else{ print ">> OK\n"; }
}

######################################################################################################

print "\n### 2) checking optional parts: \n\n";

if(!$COGSonly)
{
  # check PFAM binary and DB
  print "\n## checking optional HMMER binaries (lib/phyTools: \$ENV{'EXE_HMMPFAM'})\n";
  print "# required by get_homologues.pl -D\n";

  if($downloadOK){
    $ENV{'EXE_HMMPFAM'} = $ENV{'MARFIL'}."/bin/hmmer-3.1b2/binaries/hmmscan"
  }

  $output = `$ENV{'EXE_HMMPFAM'} 2>&1`;
  if($output =~ /Usage:/){ print ">> OK\n"; }
  else
  {
    die "<< Cannot run shipped hmmscan, please download hmmer 3.1b2 from http://hmmer.janelia.org ,\n".
      "<< install it and edit variable EXE_HMMPFAM as explained in the manual\n".
      "<<(inside set_phyTools_env in file lib/phyTools.pm) .\n".
      "<< Then re-run\n";
  }
}

if(!$noDBs || $noDBs==1)
{
  print "## checking optional PFAM library (lib/phyTools: \$ENV{'PFAMDB'})\n";
  print "# required by get_homologues.pl -D and get_homologues-est.pl -D\n";
  if(!$noDBs && ($force_unsupervised || ! -s $ENV{'PFAMDB'}.'.h3m'))
  {
    my $userword = '';
    if(!$force_unsupervised)
    {
      print "# cannot locate Pfam-A, would you like to download it now? [Y/n]\n";
      $userword = <STDIN>;
    }    

    if($force_unsupervised || $userword =~ m/Y/i)
    {
      chdir($ENV{'MARFIL'}.'/db/');
      my $ftp;
    
      if(! -s $PFAMHMMFILE)
      {
        print "# connecting to $PFAMSERVERURL ...\n";
        eval{ require Net::FTP; };
  
        if($ftp = Net::FTP->new($PFAMSERVERURL,Passive=>1,Debug =>0,Timeout=>60))
        {
          $ftp->login("anonymous",'-anonymous@') || die "# cannot login ". $ftp->message();
          $ftp->cwd($PFAMFOLDER) || warn "# cannot change working directory to $PFAMFOLDER ". $ftp->message();
          $ftp->binary();
          my $downsize = $ftp->size($PFAMHMMFILE);
          $ftp->hash(\*STDOUT,$downsize/20) if($downsize);
          printf("# downloading ftp://%s/%s/%s (%1.1fMb) ...\n",$PFAMSERVERURL,$PFAMFOLDER,$PFAMHMMFILE,$downsize/(1024*1024));
          print "# [        50%       ]\n# ";
          if(!$ftp->get($PFAMHMMFILE))
          {
            warn "# cannot download file $PFAMHMMFILE ". $ftp->message() ."\n\n";
            warn "<< You might manually download $PFAMHMMFILE from $PFAMSERVERURL/$PFAMFOLDER to any location\n".
              "<< and edit variable PFAMDB as to point to that location, as explained in the manual.\n".
              "<< Then re-run\n";
          }
        }
        else
        {
          warn "# cannot connect to $PFAMSERVERURL: $@\n\n";
          warn "<< You might manually download $PFAMHMMFILE from $PFAMSERVERURL/$PFAMFOLDER to any location\n".
            "<< and edit variable PFAMDB as to point to that location, as explained in the manual.\n".
            #"inside set_phyTools_env in file lib/phyTools.pm pointing to the destination folder.\n".
            "<< Then re-run\n";
        }
        
        print "\n";
        $ftp->quit();
      }
            
      if(-s $PFAMHMMFILE)
      {
        print "# gunzip $PFAMHMMFILE ...\n";
        system("gunzip $PFAMHMMFILE");
  
        my $hmmpress = (split(/hmmscan/,$ENV{'EXE_HMMPFAM'}))[0].'hmmpress';
        $PFAMHMMFILE =~ s/\.gz//;
        print "# pressing $PFAMHMMFILE ...\n";
        system("$hmmpress $PFAMHMMFILE");
        if(!-s $PFAMHMMFILE.'.h3m')
        {
          warn "# failed pressing $PFAMHMMFILE\n";
        }
        else
        { 
          unlink($PFAMHMMFILE);
          print ">> OK\n"; 
        }
          
        chdir($cwd);
      }
    }
    else
    {
      warn "<< You might manually download $PFAMHMMFILE from $PFAMSERVERURL/$PFAMFOLDER to any location\n".
        "<< and edit variable PFAMDB as to point to that location, as explained in the manual.\n".
        #"inside set_phyTools_env in file lib/phyTools.pm pointing to the destination folder.\n".
        "<< Then re-run\n";
    }
  }
  else{ print ">> OK\n"; }

  # check swissprot DB
  print "## checking optional SWISSPROT library (lib/phyTools: \$ENV{'BLASTXDB'})\n";
  print "# required by transcripts2cds.pl and transcripts2cdsCPP.pl\n";
  if(!$force_unsupervised && -s $ENV{'BLASTXDB'}.'.phr' && -s $ENV{'BLASTXDB'}.'.dmnd')
  {
    print ">> OK\n";
  }
  else
  {
    my $userword = '';
    if(!$force_unsupervised && $noDBs != 1)
    { 
      print "# cannot locate SWISSPROT, would you like to download it now? [Y/n]\n";
      $userword = <STDIN>; 
    }
  
    if($force_unsupervised || $noDBs==1 || $userword =~ m/Y/i)
    {
      chdir($ENV{'MARFIL'}.'/db/');
  
      my $FLATSWISSPROTFILE = $SWISSPROTFILE;
      $FLATSWISSPROTFILE =~ s/\.gz//;
     
      # remove path
      $SWISSPROTFILE = (split(/\//,$SWISSPROTFILE))[-1];
      
      if(!-s $FLATSWISSPROTFILE)
      {
        print "# connecting to $SWISSPROTSERVER ...\n";
        eval{ require Net::FTP; };
  
        if(my $ftp = Net::FTP->new($SWISSPROTSERVER,Passive=>1,Debug =>0,Timeout=>60))
        {
          $ftp->login("anonymous",'-anonymous@') || die "# cannot login ". $ftp->message();
          $ftp->cwd($SWISSPROTFOLDER) || warn "# cannot change working directory to $SWISSPROTFOLDER ". $ftp->message();
          $ftp->binary();
          my $downsize = $ftp->size($SWISSPROTFILE);
          $ftp->hash(\*STDOUT,$downsize/20) if($downsize);
          printf("# downloading ftp://%s/%s/%s (%1.1fMb) ...\n",$SWISSPROTSERVER,$SWISSPROTFOLDER,$SWISSPROTFILE,$downsize/(1024*1024));
          print "# [        50%       ]\n# ";
        
          if(!$ftp->get($SWISSPROTFILE))
          {
            warn "# cannot download file $SWISSPROTFILE ". $ftp->message() ."\n\n";
            warn "<< You might manually download $SWISSPROTFILE from $SWISSPROTSERVER/$SWISSPROTFOLDER to any location\n".
              "<< and edit variable BLASTXDB as to point to that location, as explained in the manual.\n".
              "<< Then re-run\n";
          }
  
          print "\n";
          $ftp->quit();
        }
        else
        {
          warn "# cannot connect to $SWISSPROTSERVER: $@\n\n";
          warn "<< You might manually download $SWISSPROTFILE from $SWISSPROTSERVER/$SWISSPROTFOLDER to any location\n".
            "<< and edit variable BLASTXDB as to point to that location, as explained in the manual.\n".
            #"inside set_phyTools_env in file lib/phyTools.pm pointing to the destination folder.\n".
            "<< Then re-run\n";
        }
      }      
      
      if(-s $SWISSPROTFILE)
      {
        print "# gunzip $SWISSPROTFILE ...\n"; 
        system("gunzip $SWISSPROTFILE"); 
      }
      
      if(-s $FLATSWISSPROTFILE)
      {
        executeFORMATDB_EST($FLATSWISSPROTFILE);
        executeMAKEDB($FLATSWISSPROTFILE);
               
        if(!-s $FLATSWISSPROTFILE.'.phr')
        {
          warn "# failed formatting $SWISSPROTFILE (blastx)\n";
        }
        elsif(!-s $FLATSWISSPROTFILE.'.dmnd')
        {
          warn "# failed formatting $SWISSPROTFILE (diamond)\n";
        }
        else
        {    
          system("gzip $FLATSWISSPROTFILE"); 
          print ">> OK\n"; 
        }
      }
                
      chdir($cwd);
    }
    else
    {
      warn "<< You might manually download $SWISSPROTFILE from $SWISSPROTSERVER/$SWISSPROTFOLDER to any location\n".
        "<< and edit variable BLASTXDB as to point to that location, as explained in the manual.\n".
        #"inside set_phyTools_env in file lib/phyTools.pm pointing to the destination folder.\n".
        "<< Then re-run\n";
    }
  }
}
  
# finally check whether other optional software is installed and provide installation instructions
my (@missing_packages,@missing_perl_modules);

print "## checking optional software R (lib/phyTools: \$ENV{'EXE_R'})\n";
print "# required by compare_clusters.pl, parse_pangenome_matrix.pl -s, plot_pancore_matrix.pl and\n";
print "# plot_matrix_heatmap.sh, hcluster_pangenome_matrix.sh\n";
$output = `$ENV{'EXE_R'} --version 2>&1`;
if($output =~ /R version/){ print ">> OK\n"; }
else
{
  push(@missing_packages,'R');
  print "<< missing\n";
}

print "## checking optional Perl module GD\n";
print "# required by parse_pangenome_matrix.pl -p\n";
if(eval{ require GD }){ print ">> OK\n"; }
else
{
  push(@missing_packages,'GD');
  print "<< missing\n";
}

if(@missing_packages)
{
  if(defined($SOguess))
  {
    print "\n# operating system guess: $supportedOS[$SOguess]\n";

    if($supportedOS[$SOguess] eq 'macOSX')
    {
      print "\n# Please install fink ( http://www.finkproject.org ) if required.\n";
    }

    $command = "$manager{$supportedOS[$SOguess]} install ";
    foreach my $miss (@missing_packages)
    {
      $command .= "$packages{$miss}[$SOguess] "
    }

    print "\n# Please open a terminal with root permissions (sudo, otherwise su)\n".
      "# and paste the following commands to install the missing optional parts:\n\n";
    print "$command;\n";

    if($supportedOS[$SOguess] eq 'macOSX')
    {
      print "\n# HELP while installing optional dependencies under macOSX:\n".
        "http://cran.r-project.org/bin/macosx/\nhttp://www.bugzilla.org/docs/2.16/html/osx.html\n";
    }
  }
  else
  {
    print "# cannot guess your operating system guess\n";

    foreach my $miss (@missing_packages)
    {
      $command .= "$miss, "
    }

    die "<< Sorry, cannot identify your operating system, please check the manual for instructions,\n".
      "<< install packages $command and re-run.\n";
  }
}

if(@missing_packages || @missing_perl_modules)
{
  print "\n# In case of errors please check the manual for instructions and re-run.\n";
}
else
{
  print "\n### 3) Your get_homologues kit is now fully functional\n";
}

exit(0);
