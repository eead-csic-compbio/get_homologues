package HPCluster;

# Package to manage cluster jobs from get_homologues.pl & get_homologues-est.pl

# Currently supports SGE, SLURM and LSF clusters, 
# but it should not be too dificult to add support for other systems (help welcome)

# Bruno Contreras-Moreira, Pablo Vinuesa
# 2020 CCG/UNAM, Mexico, EEAD/CSIC, Zaragoza, Spain

use strict;
require Exporter;

our @ISA = qw( Exporter );
our @EXPORT = qw( 
  read_cluster_config print_cluster_config cluster_is_available 
  submit_cluster_job check_cluster_jobs
);

# key cluster management binaries
my @CLBINS = ('SUBEXE','CHKEXE','DELEXE');

# Default SGE cluster configuration options
# can be overriden by custom config file 
my %CLUSTER_CONF;
$CLUSTER_CONF{'PATH'}   = '';     # should end with /
$CLUSTER_CONF{'TYPE'}   = 'sge';
$CLUSTER_CONF{'SUBEXE'} = 'qsub';
$CLUSTER_CONF{'CHKEXE'} = 'qstat';
$CLUSTER_CONF{'DELEXE'} = 'qdel';
$CLUSTER_CONF{'ERROR'}  = 'Eqw';  # state of failed jobs 
$CLUSTER_CONF{'QARGS'}  = '';     # queue name, resources, etc
$CLUSTER_CONF{'STIME'}  = 1;      # interval in seconds between sub commands
$CLUSTER_CONF{'CTIME'}  = 30;     # interval in seconds between stat commands

# check sample.HPC.conf for suggested LSF & slurm parameters

# Checks whether cluster config file exists and parses it.
# input: 
# 1 (string) full path to optional config file
sub read_cluster_config
{
  my ($config_file) = @_;
  if(open(CONF,"<",$config_file))
  {
    while(<CONF>)
    {
      next if(/^#/);
      chomp; 
      if(/^(\S+)\s+(\S+)/){ $CLUSTER_CONF{$1} = $2 }
    }
    close(CONF);
  }
  else 
  {
    print "# INFO: no cluster config file\n\n";
  }
}

# Prints to stdout current cluster configuration
sub print_cluster_config
{
  foreach my $conf (sort keys(%CLUSTER_CONF))
  {
    print "# $conf\t$CLUSTER_CONF{$conf}\n";
  }
  print "\n";
}

# Checks whether cluster management binaries can be used, returns 0 otherwise
# Uses system 'which'
sub cluster_is_available
{
  my ($path,$output);
  for my $bin (@CLBINS)
  {
    # concat path and binary and do system call 
    $path = `which $CLUSTER_CONF{"PATH"}$CLUSTER_CONF{$bin}`;
    chomp($path);
    if($path eq '' || $path =~ /no $CLUSTER_CONF{$bin} in/)
    {
      print "# ERROR: cannot find cluster binary $bin\n\n";
      return 0;
    } 
	else 
    {
      $output = `$CLUSTER_CONF{"PATH"}$CLUSTER_CONF{$bin} -help 2>&1`;
      if(!$output || 
        ($output !~ /usage:/i && $output !~ /use/i && $output !~ /invalid option/))
      { 
        print "# ERROR: wrong cluster binary $bin\n";
        print "$CLUSTER_CONF{'PATH'}$CLUSTER_CONF{$bin} -help\n";
        return 0; 
      }
    }
  }

  return 1;
}

# submits a cluster job, stores the assigned process id and waits STIME 
# input:
# 1 (string) job name
# 2 (string) command to be run
# 3 (string) name of output file
# 4 (string) name of work directory
# 5 reference to cluster job hash
sub submit_cluster_job
{
  my ($jobname,$command,$outfile,$dir,$ref_cluster_PIDs) = @_;

  my ($qPID,$qsubcommand) = ('','');

  if($CLUSTER_CONF{'TYPE'} eq 'lsf')
  {
    $qsubcommand = " -J n$jobname -o $outfile ";
  }
  elsif($CLUSTER_CONF{'TYPE'} eq 'slurm')
  {
    $qsubcommand = " --job-name=n$jobname -o $outfile ";
  }

  # other cluster management types could be added here with elsif

  else # default SGE
  {
    $qsubcommand = " -N n$jobname -j y -o $outfile -S /bin/bash";
  }

  $qPID = `$CLUSTER_CONF{'PATH'}$CLUSTER_CONF{'SUBEXE'} $CLUSTER_CONF{'QARGS'} $qsubcommand <<EOF
#!/bin/env bash
cd $dir
$command
EOF`;

  # parse and save process id   
  if($qPID =~ /^(\d+)\./ || 
    $qPID =~ /^Your job (\d+)/ ||
    $qPID =~ /^Job <(\d+)> is submitted/ ||
    $qPID =~ /^Submitted batch job (\d+)/ ){ $qPID = $1 } 

  # save job details associated to process id
  $ref_cluster_PIDs->{$qPID}{'command'} = $qsubcommand;
  $ref_cluster_PIDs->{$qPID}{'executable'} = $command;
  $ref_cluster_PIDs->{$qPID}{'status'} = 'sent';

  # sleep to avoid overloading
  sleep($CLUSTER_CONF{'STIME'});
}

# Checks status of cluster jobs and prints messages to stdout
# input:
# 1 (string) name of work directory
# 2 reference to cluster job hash
sub check_cluster_jobs
{
  my ($dir,$ref_PIDs) = @_;

  my ($waiting,$qPID,$newqPID,$qout) = (1);

  while($waiting)
  {
    $waiting=0;
    foreach $qPID (sort {$a<=>$b} (keys(%$ref_PIDs)))
    {
      next if($ref_PIDs->{$qPID}{'status'} eq 'deleted');

      # get status of this job
      $qout = `$CLUSTER_CONF{'PATH'}$CLUSTER_CONF{'CHKEXE'} | grep $qPID`; 
      if($qout)
      {
        if($qout =~ /\s+$CLUSTER_CONF{'ERROR'}\s+/)
        {
          # resubmit failed jobs
          $newqPID = `$CLUSTER_CONF{'PATH'}$CLUSTER_CONF{'SUBEXE'} $CLUSTER_CONF{'QARGS'} $ref_PIDs->{$qPID}{'command'} <<EOF		  
          cd $dir
          $ref_PIDs->{$qPID}{'executable'}
EOF`; 
          
          if($newqPID =~ /^(\d+)\./ || $newqPID =~ /^Your job (\d+)/){ $newqPID = $1 }
          $ref_PIDs->{$newqPID}{'command'} = $ref_PIDs->{$qPID}{'command'};
          $ref_PIDs->{$newqPID}{'executable'} = $ref_PIDs->{$qPID}{'executable'};
          $ref_PIDs->{$newqPID}{'status'} = 'sent';
          sleep($CLUSTER_CONF{'STIME'});

          # remove failed job
          system("$CLUSTER_CONF{'PATH'}$CLUSTER_CONF{'DELEXE'} $qPID");
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
      sleep($CLUSTER_CONF{'CTIME'});
    }
  }
}

1;
