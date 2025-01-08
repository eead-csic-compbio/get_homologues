# Bruno Contreras-Moreira, Pablo Vinuesa
# 2005-25 CCG-UNAM, Mexico, EEAD-CSIC, Zaragoza, Spain
# This is a collection of subroutines used in our projects,
# including primers4clades and get_homologues

package phyTools;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(
  set_phyTools_env $ERROR read_FASTA_sequence feature_is_installed check_installed_features
  extract_gene_name extract_CDSs_from_genbank extract_intergenic_from_genbank 
  extract_features_from_genbank extract_genes_from_genbank
  find_taxa_FASTA_headers find_taxa_FASTA_array_headers read_FASTA_file_array 
  same_sequence_order run_PARS check_variants_FASTA_alignment collapse_taxon_alignments
  fisher_yates_shuffle short_path factorial calc_mean calc_std_deviation 
  calc_memory_footprint warn_missing_soft
  NAME SEQ 
  );

$phyTools::ERROR = "phyError";
$phyTools::NERROR = "-999";

use strict;
use File::Basename;
use Cwd;
use File::Temp qw/ tempfile /;

# bioperl installation path if not standard
# use lib "";

eval{ require Bio::Seq } || die "phyTools: failed to import Bio::Seq, please install Bioperl (please check manual)\n@_";
eval{ require Bio::SeqIO } || die "phyTools: failed to import Bio::SeqIO, please install Bioperl (please check manual)\n@_";

#################################################################

use FindBin '$Bin';
my $MARFILPATH = "$Bin";

$phyTools::ERROR = "phyError";
$phyTools::NERROR = "-999";

# constant as indices to access read_FASTA_file_ arrays
use constant NAME => 0;
use constant SEQ  => 1;
use constant IDENTICALS => 2;
use constant NUMBER_IDENTICALS => 3;

my %feature_output = (

  # default output of binaries when corrected installed, add as required for
  # check_installed_features
  'EXE_BLASTP'=>'BLAST','EXE_BLASTN'=>'BLAST','EXE_FORMATDB'=>'database',
  'EXE_MCL'=>'usage','EXE_MULTIPARANOID'=>'args','EXE_R'=>'R is free software',
  'EXE_MAKEHASH'=>'t open output file','EXE_READBLAST'=>'Error opening','EXE_COGLSE'=>'Usage:',
  'EXE_HMMPFAM'=>'Usage:','EXE_COGTRI'=>'Error opening', #'EXE_COGTRI'=>'Runs',
  'EXE_INPARA'=>'usage','EXE_ORTHO'=>'usage','EXE_HOMOL'=>'usage','EXE_ISOFORM'=>'usage','EXE_SPLITBLAST'=>'Usage',
  'EXE_SPLITHMMPFAM'=>'Usage','EXE_MVIEW'=>'usage',
  'EXE_BLASTX_EST'=>'BLAST','EXE_TRANSDECOD_EST'=>'transcript',
  'EXE_DMNDX_EST'=>'diamond'
  );

################################################################

set_phyTools_env();

#foreach my $key (keys %ENV){print "# $key => $ENV{$key}\n";} exit;

##########################################################################################

sub set_phyTools_env
{
  if( ! defined($ENV{'MARFIL_MISSING_BINARIES'}) ) { $ENV{'MARFIL_MISSING_BINARIES'} =  ''; }
  if( ! defined($ENV{'MARFIL'}) ) { $ENV{'MARFIL'} =  $MARFILPATH .'/'; }
  
  # data 
  if( ! defined($ENV{"PFAMDB"}) ) { $ENV{"PFAMDB"} = $ENV{'MARFIL'}."db/Pfam-A.hmm"; } # HMMER3.0 not compatible with Pfam > 27
  if( ! defined($ENV{'BLASTXDB'}) ){ $ENV{'BLASTXDB'} = $ENV{'MARFIL'}."db/uniprot_sprot.fasta"; }

  # scripts
  if( ! defined($ENV{"EXE_SPLITBLAST"}) ){ $ENV{"EXE_SPLITBLAST"} = $ENV{'MARFIL'}."_split_blast.pl"; }
  if( ! defined($ENV{"EXE_SPLITHMMPFAM"}) ){ $ENV{"EXE_SPLITHMMPFAM"} = $ENV{'MARFIL'}."_split_hmmscan.pl"; }
  if( ! defined($ENV{"EXE_INPARA"}) ){ $ENV{"EXE_INPARA"} = $ENV{'MARFIL'}."_cluster_makeInparalog.pl"; }
  if( ! defined($ENV{"EXE_ORTHO"}) ){ $ENV{"EXE_ORTHO"} = $ENV{'MARFIL'}."_cluster_makeOrtholog.pl"; }
  if( ! defined($ENV{"EXE_HOMOL"}) ){ $ENV{"EXE_HOMOL"} = $ENV{'MARFIL'}."_cluster_makeHomolog.pl"; }
  if( ! defined($ENV{"EXE_ISOFORM"}) ){ $ENV{"EXE_ISOFORM"} = $ENV{'MARFIL'}."_cluster_makeIsoform.pl"; }

  # BLAST
  if( ! defined($ENV{'BLAST_PATH'}) ){ 
    $ENV{'BLAST_PATH'} = $ENV{'MARFIL'}.'bin/ncbi-blast-2.16.0+/bin/'; 
    if(!-e $ENV{'BLAST_PATH'}){ 
      $ENV{'BLAST_PATH'} = ''; # should work if in PATH
    } 
  }
  if( ! defined($ENV{'EXE_BLASTP'}) ){ $ENV{'EXE_BLASTP'} = $ENV{'BLAST_PATH'}.'blastp'; }
  if( ! defined($ENV{'EXE_BLASTN'}) ){ $ENV{'EXE_BLASTN'} = $ENV{'BLAST_PATH'}.'blastn'; }
  if( ! defined($ENV{'EXE_FORMATDB'}) ){ $ENV{'EXE_FORMATDB'} = $ENV{'BLAST_PATH'}.'makeblastdb'; }
  
  # PHYLIP PARS
  if( ! defined($ENV{'EXE_PARS'}) ){ 
    $ENV{'EXE_PARS'} = $ENV{'MARFIL'}.'bin/phylip-3.695/exe/pars'; 
    if(!-e $ENV{'EXE_PARS'}){
      $ENV{'EXE_PARS'} = 'pars'; # should work if in PATH
    }	    
  }

  # mcl
  if( ! defined($ENV{"EXE_MCL"}) ){ 
    $ENV{"EXE_MCL"} = $ENV{'MARFIL'}."/bin/mcl-14-137/src/shmcl/mcl"; 
    if(!-e $ENV{"EXE_MCL"}){
      $ENV{"EXE_MCL"} = 'mcl'; # should work if in PATH
    }	  
  }

  # HMMER
  if( ! defined($ENV{"EXE_HMMPFAM"}) ){ 
    $ENV{"EXE_HMMPFAM"} = $ENV{'MARFIL'}."/bin/hmmer-3.1b2/binaries/hmmscan --noali --acc --cut_ga "; 
    if(!-e $ENV{'MARFIL'}."/bin/hmmer-3.1b2/binaries/hmmscan"){
      $ENV{"EXE_HMMPFAM"} = 'hmmscan --noali --acc --cut_ga '; # should work if in PATH
    }	    
  }

  # COGS
  if( ! defined($ENV{"COGS_PATH"}) ){ 
    $ENV{"COGS_PATH"} = $ENV{'MARFIL'}."/bin/COGsoft";
    if(!-e $ENV{"COGS_PATH"}){
      if( ! defined($ENV{"EXE_MAKEHASH"}) ){ $ENV{"EXE_MAKEHASH"} = "COGmakehash "; }
      if( ! defined($ENV{"EXE_READBLAST"}) ){ $ENV{"EXE_READBLAST"} = "COGreadblast "; }
      if( ! defined($ENV{"EXE_COGLSE"}) ){ $ENV{"EXE_COGLSE"} = "COGlse "; }
      if( ! defined($ENV{"EXE_COGTRI"}) ){ $ENV{"EXE_COGTRI"} = "COGtriangles "; }

    } else {
      if( ! defined($ENV{"EXE_MAKEHASH"}) ){ $ENV{"EXE_MAKEHASH"} = $ENV{'COGS_PATH'}."/COGmakehash/COGmakehash "; }
      if( ! defined($ENV{"EXE_READBLAST"}) ){ $ENV{"EXE_READBLAST"} = $ENV{'COGS_PATH'}."/COGreadblast/COGreadblast "; }
      if( ! defined($ENV{"EXE_COGLSE"}) ){ $ENV{"EXE_COGLSE"} = $ENV{'COGS_PATH'}."/COGlse/COGlse "; }
      if( ! defined($ENV{"EXE_COGTRI"}) ){ $ENV{"EXE_COGTRI"} = $ENV{'COGS_PATH'}."/COGtriangles/COGtriangles "; }
    }
  }

  # MVIEW
  if( ! defined($ENV{"EXE_MVIEW"}) ){ $ENV{"EXE_MVIEW"} = $ENV{'MARFIL'}."lib/mview/bin/mview "; }
  
  # transcripts/ESTs
  if( ! defined($ENV{"EXE_TRANSDECOD_EST"}) ){ $ENV{"EXE_TRANSDECOD_EST"} = $ENV{'MARFIL'}."/lib/est/TransDecoder_r20140704/TransDecoder "; }
  if( ! defined($ENV{'BLAST_PATH_EST'}) ){ $ENV{'BLAST_PATH_EST'} = $ENV{'BLAST_PATH'}; }
  if( ! defined($ENV{'EXE_BLASTX_EST'}) ){ $ENV{'EXE_BLASTX_EST'} = $ENV{'BLAST_PATH_EST'}.'blastx'; }
  if( ! defined($ENV{'EXE_FORMATDB_EST'}) ){ $ENV{'EXE_FORMATDB_EST'} = $ENV{'BLAST_PATH_EST'}.'makeblastdb'; }
  
  # diamond 
  if( ! defined($ENV{'DMND_PATH'}) ){ 
    $ENV{'DMND_PATH'} = $ENV{'MARFIL'}.'bin/diamond-0.8.25/'; 
    if(!-e $ENV{'DMND_PATH'}){
      $ENV{'DMND_PATH'} = ''; # should work if in PATH
    }	    
  }
  if( ! defined($ENV{'EXE_DMNFT'}) ){ $ENV{'EXE_DMNFT'} = $ENV{'DMND_PATH'}.'diamond makedb'; }
  if( ! defined($ENV{'EXE_DMNDP'}) ){ $ENV{'EXE_DMNDP'} = $ENV{'DMND_PATH'}.'diamond blastp'; }
  if( ! defined($ENV{'EXE_DMNDX_EST'}) ){ $ENV{'EXE_DMNDX_EST'} = $ENV{'DMND_PATH'}.'diamond blastx'; }
  
  if( ! defined($ENV{"EXE_R"}) ){ $ENV{"EXE_R"} = "R " } # expected to be in the path
  if( ! defined($ENV{"TRAILINGZEROS"}) ){ $ENV{"TRAILINGZEROS"} = 3; } # trailing zeros for sequence identifiers, produced by read_FASTA_sequences, (3 for max 999 sequences)
}
########################################################################################

sub check_installed_features
{

  # check all needed binaries and data sources required by functions here
  # fills $ENV{"MISSING_BINARIES"} with missing binaries
  my (@to_be_checked) = @_;
  my ($check_summary,$output) = ("\nChecking required binaries and data sources, all set in phyTools.pm :\n");
  foreach my $bin (@to_be_checked)
  {
    $check_summary .= sprintf("%18s : ",$bin);
    if($ENV{$bin})
    {
      if('EXE_MAFFT_EINSI' eq $bin){ $output = `$ENV{$bin} -h 2>&1 ` }
      elsif('EXE_R' eq $bin){ $output = `$ENV{$bin} --version 2>&1 ` }
      elsif('EXE_MVIEW' eq $bin){ $output = `$ENV{$bin} -h 2>&1 ` }
      else{ $output = `$ENV{$bin} 2>&1 ` }
      if(!$output){ $output = '' }
    }
    if($output =~ /$feature_output{$bin}/)
    {
      if($bin eq 'EXE_HMMPFAM') # special check for hammer (Pfam hmm library)
      {
        if(!-s $ENV{'PFAMDB'}.'.h3m')
        {
          $output = "Pfam-A.hmm file needs to be installed/prepared (run ./install.pl)\n";
          $ENV{"MARFIL_MISSING_BINARIES"} .= "$bin,";
        }
        #elsif(!-s $ENV{'PFAMDB'}.'.h3m'){
        #  $output = "Pfam-A.hmm file should be prepared with hmmpress (same folder as $ENV{'EXE_HMMPFAM'})\n";
        #  $ENV{"MARFIL_MISSING_BINARIES"} .= "$bin,";
        #}
        else{ $output = "OK ($ENV{$bin} $ENV{'PFAMDB'})" }
      }
      else{ $output = "OK (path:$ENV{$bin})" }
    }
    else
    {
      $ENV{"MARFIL_MISSING_BINARIES"} .= "$bin,";
      $output = " needs to be installed (wrong path:$ENV{$bin})";
    }
    $check_summary .= $output."\n";
  }#die $ENV{"MARFIL_MISSING_BINARIES"};
  return $check_summary;
}

sub feature_is_installed
{

  # checks whether a passed software can be used before calling it
  my ($feature) = @_;
  my $env_missing = $ENV{"MARFIL_MISSING_BINARIES"};
  if($feature eq 'BLAST')
  {
    if($env_missing =~ /BLASTP/ || $env_missing =~ /BLASTN/ || 
      $env_missing =~ /BLASTX/ || $env_missing =~ /FORMATDB/){ return 0 }
  }
  elsif($feature eq 'DIAMOND')
  {
    if($env_missing =~ /DMNDP/ || $env_missing =~ /DMNDX/){ return 0 }
  }
  elsif($feature eq 'PFAM')
  {
    if($env_missing =~ /PFAM/){ return 0 }
  }
  elsif($feature eq 'OMCL')
  {
    if($env_missing =~ /MCL/){ return 0 }
  }
  elsif($feature eq 'PRND')
  {
    if($env_missing =~ /PARANOID/){ return 0 }
  }
  elsif($feature eq 'COGS')
  {
    if($env_missing =~ /HASH/ || $env_missing =~ /COG/ || $env_missing =~ /READBLAST/){ return 0 }
  }
  elsif($feature eq 'MAFFT')
  {
    if($env_missing =~ /MAFFT/){ return 0 }
  }
  elsif($feature eq 'R')
  {
    if($env_missing =~ /EXE_R/){ return 0 }
  }
  elsif($feature eq 'EXE_TRANSDECOD_EST')
  {
    if($env_missing =~ /HASH/ || $env_missing =~ /COG/ || $env_missing =~ /READBLAST/){ return 0 }
  }
  elsif($feature eq 'MAFFT')
  {
    if($env_missing =~ /MAFFT/){ return 0 }
  }
  elsif($feature eq 'R')
  {
    if($env_missing =~ /EXE_R/){ return 0 }
  }

  return 1;
}

# needs a pangenome_matrix.tab file as single argument to run
# Pablo May2011
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
########################################################################################

sub run_PARS
{
  # runs phylip pars from a binary matrix file and returns a tree and a log generated by PHYLIP PARS
  # optionally takes a hash of label equivalences to rename tree leaves
  # WARNING: in/out files are always named the same, so cannot run more than one process simultaneosuly
  # in the same directory!
  my($infile,$treefile,$treelogfile,$ref_hash_labels) = @_;
  my ($tmp_infile,$tmp_outfile,$tmp_outtree)  = ('infile','outfile','outtree');
  my ($treeline,$parsOK) = ('',0);

  # 1) copy input file to tmp input
  system("cp $infile $tmp_infile");
  my($paramsfh,$paramsfilename) = tempfile();
  print $paramsfh "j\n7\n50\ny\n";
  close($paramsfh);

  # 2) run pars
  open(PARS,"$ENV{'EXE_PARS'} < $paramsfilename |") || die "# run_PARS: cannot run $ENV{'EXE_PARS'}\n";
  while(<PARS>)
  {

    #print;
    if(/^Output written to file/){ $parsOK++ }
    elsif(/^Tree also written onto file/){ $parsOK++ }
  }
  close(PARS);

  # 3) create log file
  rename($tmp_outfile,$treelogfile);

  # 4) produce tree file
  if($ref_hash_labels)
  {

    # parse $tmp_outtree get rid of \n and \s, and write the single treeline to file
    open(TMP,$tmp_outtree) || die "run_PARS cannot open outfile $tmp_outtree: $!\n";
    while(<TMP>)
    {
      $_ =~ s/\n//g;
      $treeline .= $_;
    }
    close(TMP);

    # create labelled tree
    my $labelled_newick = add_labels2newick_tree( $treeline, $ref_hash_labels );
    open(LABELTREE,">$treefile") || die "# run_PARS : cannot create $treefile\n";
    print LABELTREE $labelled_newick;
    close(LABELTREE);
    unlink($tmp_outtree);
  }
  else{ rename($tmp_outtree,$treefile); }

  # 5) clean temporary files
  unlink($tmp_infile,$tmp_outtree);
}

# edited in Dec2014 to support compressed files
sub read_FASTA_sequence
{
  # if $remove_gap_col assumes FASTA file is actually a multiple alignment
  # in FASTA format
  # if $skipbadCDSs ignores nt sequences with length%3 > 0
  my ( $infile, $remove_gap_cols, $skipbadCDSs, $skipidentical ) = @_;
  my (%FASTA,$magic,$name,$seq,$n_of_sequences,$length,$maxlength,$pos,$isgap,$seqid);
  $n_of_sequences = $maxlength = 0;

  # check input file format and open it accordingly
  open(INFILE,$infile) || die "# read_FASTA_sequence: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    if(!open(FASTA,"gzip -dc $infile |"))
    {
      die "# read_FASTA_sequence: cannot read GZIP compressed $infile $!\n"
        ."# please check gzip is installed\n";
    }
  }
  elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    if(!open(FASTA,"bzip2 -dc $infile |"))
    {
      die "# read_FASTA_sequence: cannot read BZIP2 compressed $infile $!\n"
        ."# please check bzip2 is installed\n";
    }
  }
  else{ open(FASTA,"<$infile") || die "# read_FASTA_sequence: cannot read $infile $!\n"; }

  while(<FASTA>)
  {
    if(/^\>/)
    {
      $name = $_;
      $n_of_sequences++;
      $seqid = sprintf("%0$ENV{'TRAILINGZEROS'}d",$n_of_sequences); #print "#$seqid#\n";
      $FASTA{$seqid}{'NAME'} = $name;
    }
    else
    {
      $FASTA{$seqid}{'SEQ'} .= $_;
      $FASTA{$seqid}{'SEQ'} =~ s/[\s|\n]//g;
      $length = length($FASTA{$seqid}{'SEQ'});
      if($length > $maxlength){ $maxlength = $length; }
    }
  }
  close(FASTA);
  if($skipbadCDSs)
  {
    # check length is divisible by 3 and that no inframe STOP codons exist
    my %goodCDSs_FASTA;
    foreach $seq (keys(%FASTA))
    {
      my $ntseq = Bio::Seq->new( -display_id => 'tmp', -seq => $FASTA{$seq}{'SEQ'} );
      my $protseq = $ntseq->translate()->seq(); chop $protseq;
      if($protseq =~ /\*/)
      {
        print "# read_FASTA_sequence : skipped CDS sequence (inframe STOP codon) $FASTA{$seq}{'NAME'}\n";
        next;
      }
      if($FASTA{$seq}{'SEQ'} !~ /[QEIPDFHLV]/ && length($FASTA{$seq}{'SEQ'})%3)
      {
        my $choppednts = 0;
        while(length($FASTA{$seq}{'SEQ'})%3){ chop $FASTA{$seq}{'SEQ'}; $choppednts++; }
        print "# read_FASTA_sequence : chopped CDS sequence (3' $choppednts nts) $FASTA{$seq}{'NAME'}\n";
      }

      #if(length($protseq) < 100){ print "mira: $protseq\n"; }
      $goodCDSs_FASTA{$seq} = $FASTA{$seq};
    }
    %FASTA = %goodCDSs_FASTA;
  }
  if($skipidentical)
  {
    my (%nrFASTA,%identical,$seq2);
    foreach $seq (keys(%FASTA))
    {
      $FASTA{$seq}{'IDENTICALS'} = '';
      next if($identical{$seq});
      foreach $seq2 (keys(%FASTA))
      {
        next if($seq == $seq2);
        if($FASTA{$seq}{'SEQ'} eq $FASTA{$seq2}{'SEQ'} ||        # same length
          $FASTA{$seq}{'SEQ'} =~ /$FASTA{$seq2}{'SEQ'}/ || # different length
          $FASTA{$seq2}{'SEQ'} =~ /$FASTA{$seq}{'SEQ'}/)
        {
          $identical{$seq2} = $seq;
          $FASTA{$seq}{'IDENTICALS'} .= "$seq2,"; #print "mira $seq $seq2\n";
          $FASTA{$seq}{'NUMBER_IDENTICALS'}++;
        }
      }
    }
    foreach $seq (keys(%FASTA))
    {
      if($identical{$seq})
      {
        print "# read_FASTA_sequence : skipped sequence identical to $identical{$seq}: $FASTA{$seq}{'NAME'}\n";
        next;
      }
      $nrFASTA{$seq} = $FASTA{$seq};

      # keep track of identical sequences
      if($FASTA{$seq}{'NUMBER_IDENTICALS'})
      {
        chomp($nrFASTA{$seq}{'NAME'});
        $nrFASTA{$seq}{'NAME'} .= " | identical sequences=$FASTA{$seq}{'NUMBER_IDENTICALS'}\n";
      }
    }
    %FASTA = %nrFASTA;
  }

  #if($remove_gap_cols && $remove_gap_cols == 1)
  if($remove_gap_cols)
  {
    for($pos=0;$pos < $length; $pos++)
    {
      $isgap=0;
      foreach $seq (keys(%FASTA))
      {
        if(substr($FASTA{$seq}{'SEQ'},$pos,1))
        {
          if(substr($FASTA{$seq}{'SEQ'},$pos,1) eq '-'){ $isgap=1; last; }
        }
        else{ $isgap=1; last }
      }
      if($isgap == 0)
      {
        foreach $seq (keys(%FASTA))
        {
          $FASTA{$seq}{'NOGAPSEQ'} .= substr($FASTA{$seq}{'SEQ'},$pos,1);
        }
      }
    }
    foreach $seq (keys(%FASTA))
    {
      $FASTA{$seq}{'SEQ'} = $FASTA{$seq}{'NOGAPSEQ'}
    }
  }
  return %FASTA;
}

# edited in Dec2014 to support compressed files
# added in 2011 as arrays take less RAM and are faster to browse:
# http://bioinfoperl.blogspot.com/2011/02/estructuras-de-datos-basicas-en-perl.html
sub read_FASTA_file_array
{
  # in FASTA format
  # returns a reference to a 2D array for 4 secondary keys: NAME,SEQ,IDENTICALS,NUMBER_IDENTICALS
  # first valid index (first sequence) is '0'
  my ( $infile, $skipamino, $skipidentical ) = @_;
  my (@FASTA,$name,$seq,$n_of_sequences,$magic);
  $n_of_sequences = -1;

  # check input file format and open it accordingly
  open(INFILE,$infile) || die "# read_FASTA_sequence_array: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    if(!open(FASTA,"gzip -dc $infile |"))
    {
      die "# read_FASTA_sequence_array: cannot read GZIP compressed $infile $!\n"
        ."# please check gzip is installed\n";
    }
  }
  elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    if(!open(FASTA,"bzip2 -dc $infile |"))
    {
      die "# read_FASTA_sequence_array: cannot read BZIP2 compressed $infile $!\n"
        ."# please check bzip2 is installed\n";
    }
  }
  else{ open(FASTA,"<$infile") || die "# read_FASTA_sequence_array: cannot read $infile $!\n"; }

  while(<FASTA>)
  {
    next if(/^$/ || /^#/);
    if(/^\>(.*?)[\n\r]/)
    {
      $n_of_sequences++; # first sequence ID is 0
      $name = $1; #print "# $n_of_sequences $name\n";
      $FASTA[$n_of_sequences][NAME] = $name;
    }
    elsif($n_of_sequences>-1)
    {
      $_ =~ s/[\s|\n|\-|\.]//g;
      $FASTA[$n_of_sequences][SEQ] .= $_; # conserve case
    }
  }
  close(FASTA);

  if($n_of_sequences == -1){ return \@FASTA }
  if($skipamino || $skipidentical)
  {
    my ($n_of_nr_sequences,@nrFASTA,%identical) = 0;
    if($skipidentical)
    {
      foreach $seq ( 0 .. $n_of_sequences-1 )
      {
        $FASTA[$seq][IDENTICALS] = '';
        next if($identical{$seq});
        foreach my $seq2 ( $seq + 1 .. $n_of_sequences )
        {
          if($FASTA[$seq][SEQ] eq $FASTA[$seq2][SEQ] ||  # same length
            $FASTA[$seq][SEQ] =~ /$FASTA[$seq2][SEQ]/ || # different length
            $FASTA[$seq2][SEQ] =~ /$FASTA[$seq][SEQ]/)
          {
            $identical{$seq2} = $seq;
            $FASTA[$seq][IDENTICALS] .= "$seq2,"; #print "mira $seq $FASTA[$seq][IDENTICALS]\n";
            $FASTA[$seq][NUMBER_IDENTICALS]++;
          }
        }
      }
    }
    foreach $seq ( 0 .. $n_of_sequences )
    {
      # skip protein sequences
      if($FASTA[$seq][SEQ] =~ /[QEIPDFHLV]/i)
      {
        print "# read_FASTA_sequence_array : skipped amino acid sequence ($FASTA[$seq][NAME])\n";
        next;
      }
      if($identical{$seq})
      {
        print "# read_FASTA_sequence_array : skipped identical sequence: ".
          "$FASTA[$seq][NAME]\n";
        next;
      }
      $nrFASTA[$n_of_nr_sequences] = $FASTA[$seq];

      # keep track of identical sequences
      if($FASTA[$seq][NUMBER_IDENTICALS])
      {
        chomp($nrFASTA[$n_of_nr_sequences][NAME]);
        $nrFASTA[$n_of_nr_sequences][NAME] .=
          " | identical sequences=$nrFASTA[$seq][NUMBER_IDENTICALS]";
      }
      $n_of_nr_sequences++;
    }
    return \@nrFASTA;
  }
  else{ return \@FASTA }
}

sub same_sequence_order
{
  # Takes two refs to arrays of same length made by read_FASTA_sequence_array,
  # 1st nucleotide CDS and 2nd peptide CDS.
  # Returns 1 if sequences are in same order, otherwise 0
  # Actually translates sequences, does not use sequence names. 
  # Note: nucleotide CDS with N residues are not compared 

  my ($refCDSnt, $refCDSaa) = @_;

  my ($seq, $translated, $parsed);

  foreach $seq (0 .. $#{$refCDSnt})
  {
    next if($refCDSnt->[$seq][SEQ] =~ /N/i);

    $translated = Bio::Seq->new(-seq => $refCDSnt->[$seq][SEQ])->translate()->seq();
    $parsed = $refCDSaa->[$seq][SEQ];
    $translated =~ s/\*$//g;
    $parsed =~ s/\*$//gi;
   
    if($translated =~ m/\*/) {
      print "# same_sequence_order : translation of $refCDSnt->[$seq][NAME] contains internal STOP codons (*): ".
        "$translated\n";
      return 0;
    }   
   
    if($parsed ne $translated && 
      $parsed !~ /\Q$translated\E/ && $translated !~ /\Q$parsed\E/) {
      print "# same_sequence_order : sequence #$seq does not match: ".
        "$refCDSnt->[$seq][SEQ]\n$translated\n$refCDSaa->[$seq][SEQ]\n";
      return 0;
    }
  }

  return 1;
}

sub find_taxa_FASTA_headers
{
  # takes a read_FASTA_sequence-returned hash ref and finds all present taxa, parsed as [strings]
  # expects [taxon name] to the righttmost side, and applies some naive filters
  # returns a hash with taxa names as keys and occurences as values
  
  my($ref_FASTA) = @_;
  my ($seq,$taxon,$header,%TAXA);
  foreach $seq (keys(%{$ref_FASTA}))
  {
    # assume last [] is the taxon name, to avoid problems with headers such as
    # [acyl-carrier-protein] synthase I [Buchnera aphidicola str. Bp (Baizongia pistaciae)]
    $header = reverse($ref_FASTA->{$seq}{'NAME'});
    while($header =~ /(\].+?\[)/g)
    {
      $taxon = reverse($1);

      # adhoc rules to filter out common annotations such as:
      # [5'-phosphate] , [3-hydroxymyristoyl] ,
      # [acyl-carrier-protein] , [NAD(P)H] , [Cu-Zn]
      next if(substr($taxon,1,1) !~ /[A-Z]/ || $taxon !~ /\w+\s+\w+/);
      $TAXA{$taxon}{'SIZE'}++;
      push(@{$TAXA{$taxon}{'MEMBERS'}},$seq);
      last;
    }
  }
  return %TAXA;
}

sub find_taxa_FASTA_array_headers
{
  # takes a read_FASTA_file_array-returned array ref and finds all present taxa, parsed as [strings]
  # expects [taxon name] to the rightmost side, and applies some naive filters
  # returns a hash with taxa names as keys and SIZE and MEMBERS as (hash) values
  # if $single_word_taxons_are_OK is passed, [labels] are parsed as valid taxon names, as opposed to
  # default [Genus species ...]
  
  my($ref_FASTA,$single_word_taxons_are_OK) = @_;
  my ($seq,$taxon,$header,%TAXA);
  
  foreach $seq (0 .. $#{ $ref_FASTA } )
  {
    # assume last [] is the taxon name, to avoid problems with headers such as
    # [acyl-carrier-protein] synthase I [Buchnera aphidicola str. Bp (Baizongia pistaciae)]
    $header = reverse($ref_FASTA->[$seq][NAME]);
    while($header =~ /(\].+?\[)/g)
    {
      $taxon = reverse($1);

      # adhoc rules to filter out common annotations such as:
      # [5'-phosphate] , [3-hydroxymyristoyl] ,
      # [acyl-carrier-protein] , [NAD(P)H] , [Cu-Zn]
      #[Bdistachyonv1.ABR2.1.primaryTrs.cds.fna.gz]

      next if(!$single_word_taxons_are_OK &&
        (substr($taxon,1,1) !~ /[A-Z]/ || $taxon !~ /\w+\s+\w+/)); #print ">> $taxon <<\n";
          
      $TAXA{$taxon}{'SIZE'}++; 
      push(@{$TAXA{$taxon}{'MEMBERS'}},$seq);
      last;
    }
  }
  return %TAXA;
}

sub add_labels2newick_tree
{
  # assumes a hash f hashes with {001}{NAME} = "ftz Ecoli"
  #(0000000002:20.50,((0000000011:8.67,(((0000000008:14.33,0000000007:12.67):4.00,0000000005:12.33):5.33,(0000000009:2.00,(0000000010:5.00,0000000006:7.00,(0000000004:14.00,0000000003:0.00):1.00):1.00):1.33):2.00):7.83,0000000001:43.50):6.00,0000000000:34.50)[0.2000];
  #(((0000000011:9.67,(((0000000008:14.83,0000000007:12.17):4.00,0000000005:11.83):4.83,(0000000006:6.00,(0000000010:5.00,0000000009:3.00,(0000000004:14.00,0000000003:0.00):1.00):1.00):2.33):2.00):7.33,0000000002:20.00):6.00,0000000001:44.00,0000000000:34.00)[0.2000];
  #(((0000000011:9.67,(((0000000008:14.83,0000000007:12.17):4.00,0000000005:11.83):4.83,(0000000006:6.50,(0000000009:2.50,(0000000010:5.00,(0000000004:14.00,0000000003:0.00):1.00):0.50):1.00):1.83):2.00):7.33,0000000002:20.00):6.00,0000000001:44.00,0000000000:34.00)[0.2000];
  #(((0000000011:8.87,(((0000000008:14.43,0000000007:12.57):4.00,0000000005:12.23):5.23,(0000000009:2.00,(0000000010:5.00,0000000006:7.00,(0000000004:14.00,0000000003:0.00):1.00):1.00):1.53):2.00):7.73,0000000002:20.40):6.00,0000000001:43.60,0000000000:34.40)[0.2000];
  #(((0000000011:9.00,(((0000000008:14.50,0000000007:12.50):4.00,0000000005:12.17):5.17,(0000000009:2.17,(0000000006:6.83,(0000000010:5.00,(0000000004:14.00,0000000003:0.00):1.00):0.17):1.00):1.50):2.00):7.67,0000000002:20.33):6.00,0000000001:43.67,0000000000:34.33)[0.2000];

  my ($fully_labelled_tree,$ref_label_hash) = @_;
  my ($label,$oldlabel,$taxon);
  foreach my $seq (keys(%{$ref_label_hash}))
  {
    $taxon = '';
    $label = $ref_label_hash->{$seq}{'NAME'};
    $label =~ s/^>//g;
    $label = (split(/#\d+\|/,$label))[0]; # ignora primers  #164|352|0|1|1|0|55.4|57.8|AGTCCAAG...
    $label =~ s/<<.*?>>//g;
    $label =~ s/\s+|;|\-|,|:|\*|\.|\(|\)|#/_/g;
    $label =~ s/_{1,}/_/g;
    $label =~ s/==\S*?==//g;
    $label =~ s/$seq\_//g; #print "mira:$label\n";
    if($label =~ /(\[.*?\])/){ $taxon = $1 }
    if($taxon && substr($label,0,100) !~ /$taxon/){ $label = substr($label,0,15) . '...' . $taxon; } #5May09 PV ~ if($taxon && substr($label,0,25) !~ /$taxon/){ $label = substr($label,0,15) . '...' . $taxon; }
    else{ $label = substr($label,0,100) . '...'; } #print  "mira:$label\n";

    if($fully_labelled_tree =~ /($seq\__\d+)/)
    {
      $oldlabel = $1;

      #$fully_labelled_tree =~ s/$oldlabel/$oldlabel $label/;
      $fully_labelled_tree =~ s/$oldlabel/$label/g;
    }
    elsif($fully_labelled_tree =~ /($seq):/)
    {
      $oldlabel = $1;

      #$fully_labelled_tree =~ s/$oldlabel:/$oldlabel $label:/g;
      $fully_labelled_tree =~ s/$oldlabel:/$label:/g;
    }
  }

  #die $fully_labelled_tree;
  return join(";\n",split(/;/,$fully_labelled_tree));
}

sub extract_intergenic_from_genbank
{
  # takes a genbank input file and creates a FNA file containing all intergenic sequences found
  # inspired by http://bioperl.org/pipermail/bioperl-l/2006-March/021065.html
  # use $min_intergenic_size = 0 to skip minimal length tests
  # use $max_intergenic_size = 0 to skip maximum length test
  # use $length_flanking_ORFs > 0 if you want to cut oligonucleotides from both flanking ORFs to be used as PCR anchors
  
  my ($infile,$out_intergenic_file,$min_intergenic_size,$max_intergenic_size,$length_flanking_ORFs) = @_;
  if(!defined($length_flanking_ORFs) || $length_flanking_ORFs < 0){ $length_flanking_ORFs = 0 }
  my ($n_of_intergenic,$gi,$start,$end,$length,$strand,$genename,$taxon,$magic,$in) = (0);
  
  # check input file and open it accordingly
  open(INFILE,$infile) || die "# read_FASTA_sequence_array: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    $in = new Bio::SeqIO(-file => "gzip -dc $infile |", -format => 'genbank' )
  }
  elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    $in = new Bio::SeqIO(-file => "bzip2 -dc $infile |", -format => 'genbank' )
  }
  else{ $in = new Bio::SeqIO(-file => $infile, -format => 'genbank' ) }

  # first create tmp outfile, in case sub does not finish cleanly
  open(FNA,">$out_intergenic_file.tmp") || die "# extract_intergenic_from_genbank : cannot create $out_intergenic_file.tmp\n";
  while( my $seq = $in->next_seq) # scan all sequences found in $input
  {
    $seq->alphabet('dna');
    my ($gbaccession,$sequence,$gen,@genes) = ( $seq->accession() );
    if($gbaccession eq 'unknown'){ $gbaccession = $seq->display_id() } # prokka-compatible
 
    $sequence = $seq->primary_seq()->seq() || 'empty, need a full genbank entry!';
    $taxon = '';
    for my $f ($seq->get_SeqFeatures)
    {
      if($f->primary_tag =~ /CDS|rRNA|tRNA/) # ~genes
      {
        #GI deprecated on Sept2016 https://www.ncbi.nlm.nih.gov/news/03-02-2016-phase-out-of-GI-numbers/
        $gi = $genename = ''; # compatible con subrutina extract_CDS_from_genbank
        if($f->has_tag('db_xref'))
        {
          my $crossrefs = join(',',sort $f->each_tag_value('db_xref'));
          #if($crossrefs =~ /(GI\:\d+)/){ $gi = $1 }
          if($crossrefs =~ /SEED:fig\|([\w\.]+)/){ $gi = $1; } # Rast support Jul2016
        }
        
        if($f->has_tag('protein_id'))
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('protein_id')); #print "$gi\n";
        }

        if($gi eq '' && $f->has_tag('locus_tag'))
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('locus_tag'));
        }

        if($f->has_tag('gene'))
        {
          $genename = join(',',sort $f->each_tag_value('gene'));
        }
        else
        {
          if($f->has_tag('product'))
          {
            $genename = join(',',sort $f->each_tag_value('product'));
            if(length($genename)>20){ $genename = substr($genename,0,20).'..' } # avoid long product names
          }
        }
        next if($gi eq '');
        $start = $f->start(); # 1..n naturals
        $end   = $f->end();
        $strand = $f->location()->strand();
        push(@genes,[$start,$end,$gi,$strand,$genename]);
      }
      elsif($f->primary_tag() =~ /source/)
      {
        if($f->has_tag('organism'))
        {
          foreach my $element ($f->each_tag_value('organism')){ $taxon .= "[$element],"; }
          chop $taxon;
        }
      }
    }

    # calcular ORFs vecinos y corta las secuencias intergenicas
    for($gen=1;$gen<scalar(@genes);$gen++)
    {
      next if($genes[$gen-1][1] >= $genes[$gen][0]); # evita solapar, asume genes estan en orden en GenBank
      next if( defined($min_intergenic_size) && $min_intergenic_size > 0 &&
        $genes[$gen][0]-$genes[$gen-1][1]+1 < $min_intergenic_size); # skip short intergenes
      next if( defined($max_intergenic_size) && $max_intergenic_size > 0 &&
        $genes[$gen][0]-$genes[$gen-1][1]+1 > $max_intergenic_size); # skip long intergenes
      next if($genes[$gen-1][1]-$genes[$gen-1][0]+1<$length_flanking_ORFs); # skip if neighbor ORFs are too short
      next if($genes[$gen][1]-$genes[$gen][0]+1<$length_flanking_ORFs);
      $n_of_intergenic++;

      # calculate intergene boundaries for header (in 1 .. n naturals)
      $start  = $genes[$gen-1][1] + 1 - $length_flanking_ORFs; # first intergene nt - flank
      $end    = $genes[$gen][0] - 1 + $length_flanking_ORFs; # last intergene nt + flank
      $length = $end - $start + 1;
      print FNA ">intergenic$n_of_intergenic|$gbaccession|coords:$start..$end|length:$length|".
        "neighbours:$genes[$gen-1][2]($genes[$gen-1][3]),$genes[$gen][2]($genes[$gen][3])|".
        "neighbour_genes:$genes[$gen-1][4],$genes[$gen][4]|".
		  basename($infile)."|$taxon\n";

      # calculate sequence boundaries for substr, which starts in 0
      #print FNA lc(substr($sequence,$start-1,$length))."\n"; # in one go!
      # left ORF upper case
      print FNA substr($sequence,$genes[$gen-1][1]-$length_flanking_ORFs,$length_flanking_ORFs);

      # intergene lower case
      print FNA lc(substr($sequence,$genes[$gen-1][1],$genes[$gen][0]-$genes[$gen-1][1]-1));

      # right ORF upper case
      print FNA substr($sequence,$genes[$gen][0]-1,$length_flanking_ORFs)."\n";
    }
  }
  close(FNA);

  # make final outfiles only when process finishes cleanly
  rename( "$out_intergenic_file.tmp" , $out_intergenic_file );
  return $n_of_intergenic;
}

sub extract_features_from_genbank
{
  # takes a genbank input file and creates a single FASTA nucleotide file containing all features
  # defined in argument $features_to_parse or $ENV{'GBKFEATURES'}
  # if $out_feature_file == 0 does not create any outfile
  # returns reference to features parsed
  
  #     gene            230511..230630
  #                     /locus_tag="ECSE_5SrRNA01"
  #     rRNA            230511..230630
  #                     /locus_tag="ECSE_5SrRNA01"
  #                     /product="5S ribosomal RNA"
  #     gene            230683..230759
  #                     /locus_tag="ECSE_tRNA03"
  #     tRNA            230683..230759
  #                     /locus_tag="ECSE_tRNA03"
  #                     /product="tRNA-Asp"
  #     gene            230922..231725
  #                     /locus_tag="ECSE_0203"
  #     CDS             230922..231725
  #                     /locus_tag="ECSE_0203"
  #                     /product="2,5-diketo-D-gluconate reductase"
  
  my ($infile,$out_feature_file,$features_to_parse) = @_;
  if(!$features_to_parse){ $features_to_parse = $ENV{'GBKFEATURES'} }
  my ($gi,$genename,$coords,$gbaccession,$header,$genelength,%already_seen);
  my ($start,$end,$strand,$source,$featseq,$rev,$magic,$in);
  my ($taxon,$strain)=('','');
  
  # check input file and open it accordingly
  open(INFILE,$infile) || die "# read_FASTA_sequence_array: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    $in = new Bio::SeqIO(-file => "gzip -dc $infile |", -format => 'genbank' )
  }
  elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    $in = new Bio::SeqIO(-file => "bzip2 -dc $infile |", -format => 'genbank' )
  }
  else{ $in = new Bio::SeqIO(-file => $infile, -format => 'genbank' ) }

  # first create tmp outfile, in case sub does not finish cleanly
  if($out_feature_file)
  {
    open(FNA,">$out_feature_file.tmp") ||
      die "# extract_features_from_genbank : cannot create $out_feature_file.tmp\n";
  }
  while( my $seq = $in->next_seq) # scan all sequences found in $input
  { 
    $seq->alphabet('dna'); 
    $gbaccession = $seq->accession();
    if($gbaccession eq 'unknown'){ $gbaccession = $seq->display_id() } # prokka-compatible
    
    $taxon = $coords = $genelength = $source = '';
    for my $f ($seq->get_SeqFeatures)
    {
      if($f->primary_tag() =~ /source/)
      {
        $source = $f->end();
        if($f->has_tag('organism'))
        {
          foreach my $element ($f->each_tag_value('organism')){ $taxon .= "[$element],"; }
          chop $taxon;
        }
      }
      elsif($strain eq '' && $f->has_tag('strain'))
      {
        $strain = join(',',sort $f->each_tag_value('strain'));
      }
      elsif($features_to_parse =~ $f->primary_tag())
      {
        $gi = $genename = $featseq = $rev = ''; # compatible con subrutina extract_CDS_from_genbank
        $coords = $f->spliced_seq();

        if($f->location->isa('Bio::Location::SplitLocationI'))
        {
          ## cut DNA sequence of this feature by appending exons in the proper orientation
          my ($exonseq,$splitloc,%strands);
          for $splitloc ($f->location->sub_Location(0))
          {
            $strand = $splitloc->strand();
            $strands{$strand}++;
          }

          # two splicing scenarios:
          # A: cut and append exons first, finally take rc if requested
          if(scalar(keys(%strands)) == 1)
          {
            for $splitloc ($f->location->sub_Location(0))
            {
              $exonseq = $seq->subseq($splitloc->start,$splitloc->end);
              $featseq .= $exonseq;
            }

            if($strand == -1)
            {
              $featseq =~ tr/acgtnACGTN/tgcanTGCAN/;
              $featseq = reverse($featseq);
            }
          }
          else # B: cut and take rc of each exon and then append
          {
            for $splitloc ($f->location->sub_Location(0))
            {
              $exonseq = $seq->subseq($splitloc->start,$splitloc->end);
              if($splitloc->strand == -1)
              {
                $exonseq =~ tr/acgtnACGTN/tgcanTGCAN/;
                $exonseq = reverse($exonseq);
              }
              $featseq .= $exonseq;
            }
          }
        }
        else
        {
          $featseq = $coords->{'seq'};
        }

        if($gi eq '' && $f->has_tag('locus_tag'))
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('locus_tag'));
        }
        elsif($gi eq '' && $f->has_tag('db_xref'))
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('db_xref'));
        }

        if($f->has_tag('gene'))
        {
          $genename = join(',',sort $f->each_tag_value('gene'));
          $genename =~ s/\s+/_/g;
        }
        elsif($f->has_tag('product'))
        {
          $genename = join(',',sort $f->each_tag_value('product'));
          if(length($genename)>20){ $genename = substr($genename,0,20).'..' } # avoid long product names
          $genename =~ s/\s+/_/g;
        }

        next if(!$coords->{'seq'});
        
        $start = $f->start(); # 1..n naturals
        $end   = $f->end();
        $strand = $f->location()->strand();
        if(!$gi) # in cases such as nonCDS items in Smelilloti Jul2012
        {
          $gi = "ID:$gbaccession:$start-$end";
        }
        if(!$already_seen{$source}{"$gi,$start,$end,$strand"})
        {
          if($out_feature_file)
          {
            $genelength=length($featseq);
            $header = $gi." |".$taxon."|".$strain."|".basename($infile)."|".$genename."|".$genelength."|".
              "$gbaccession($source):$start-$end:$strand ".$coords->desc();
            print FNA ">$header\n$featseq\n" if($out_feature_file);
          }
          $already_seen{$source}{"$gi,$start,$end,$strand"} = 1;
        }
      }
    }
  }

  if($out_feature_file)
  {
    close(FNA);

    # make final outfiles only when process finishes cleanly
    rename( "$out_feature_file.tmp" , $out_feature_file );
  }
  return \%already_seen;
}

sub extract_CDSs_from_genbank
{
  # takes a genbank input file and creates two FASTA files containing all CDSs in
  # aminoacid and dna sequences, respectively
  # returns number of CDSs found

  my ($infile,$out_prot_file,$out_dna_file,$verbose) = @_;
  my ($gene,$gi,$crossrefs,$gbaccession,$header,$genelength,$CDScoords,$CDSseq,$rev);
  my ($start,$end,$strand,$source,$protsequence,$magic,$in,%already_seen,%genes);
  my $n_of_CDS = 0;
  my ($taxon,$strain)=('','');
  
  if($verbose){ print "# extract_CDSs_from_genbank : extracting file $infile...\n" }
  
  # check input file and open it accordingly
  open(INFILE,$infile) || die "# read_FASTA_sequence_array: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    $in = new Bio::SeqIO(-file => "gzip -dc $infile |", -format => 'genbank' )
  }
  elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    $in = new Bio::SeqIO(-file => "bzip2 -dc $infile |", -format => 'genbank' )
  }
  else{ $in = new Bio::SeqIO(-file => $infile, -format => 'genbank' ) }
  
  # parse all CDS sequences and store their neighbors
  while( my $seq = $in->next_seq) # scan all sequences found in $input
  {
    # make sure it is handled as DNA
    #http://bioperl-l.bioperl.narkive.com/V9J9amXL/parseing-trna-sequences-out-of-genbank-file
    $seq->alphabet('dna'); 
    
    $gbaccession = $seq->accession();
    if($gbaccession eq 'unknown'){ $gbaccession = $seq->display_id() } # prokka-compatible
    
    $source = '';
    foreach my $f ($seq->get_SeqFeatures)
    {
      if($f->primary_tag() =~ /source/)
      {
        $source = $f->end();
        if($taxon eq '' && $f->has_tag('organism'))
        {
          foreach my $element ($f->each_tag_value('organism')){ $taxon .= "[$element],"; }
          chop $taxon;
        }
        if($strain eq '' && $f->has_tag('strain'))
        {
          $strain = join(',',sort $f->each_tag_value('strain'));
        }
      }
      elsif($f->primary_tag() =~ /CDS/)
      {
        next if $f->has_tag('pseudo'); # skip pseudogenes altogether

        #GI deprecated on Sept2016 https://www.ncbi.nlm.nih.gov/news/03-02-2016-phase-out-of-GI-numbers/
        $gene=$gi=$crossrefs=$genelength=$protsequence=$CDSseq=$rev=''; 
        $CDScoords = $f->spliced_seq();

        if($f->location->isa('Bio::Location::SplitLocationI'))
        {
          #LOCUS       NC_002505            2961149 bp    DNA     circular BCT 24-MAY-2010
          #ACCESSION   NC_002505
          #VERSION     NC_002505.1  GI:1564003
          #CDS      complement(join(552383..553567,553567..554145))

          #print $f->start() .'..'.$f->end()."\n";
          for my $location ( $f->location->sub_Location() )
          {

            #print $location->start.'..'.$location->end.' '.$location->strand."\n";
            $CDSseq .= $seq->subseq($location->start,$location->end);
            if($location->strand == -1){ $rev = 1 }
          }

          if($rev)
          {
            $CDSseq =~ tr/acgtnACGTN/tgcanTGCAN/;
            $CDSseq = reverse($CDSseq);
          } #print "$CDSseq\n";
        }
        else
        {
          $CDSseq = $CDScoords->{'seq'};
        }
        
        if($f->has_tag('db_xref'))
        {
          $crossrefs = join(',',sort $f->each_tag_value('db_xref'));
          #if($crossrefs =~ /(GI\:\d+)/){ $gi = $1; $crossrefs =~ s/$gi// }
          if($crossrefs =~ /SEED:fig\|([\w\.]+)/ || 
            $crossrefs =~ /RAST\d*:fig\|([\w\.]+)/){ $gi = $1; } # RAST support
          next if($crossrefs =~ /PSEUDO:/); # no sabemos si es universal, funciona con Bradyrizobium_ORS278.gb
        }
        
        if($f->has_tag('protein_id'))
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('protein_id')); #print "$gi\n";
        }        
        if($f->has_tag('translation'))
        {
          $protsequence = join('',sort $f->each_tag_value('translation'));
        }
        if($gi eq '' && $f->has_tag('locus_tag'))
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('locus_tag')); #print "$gi\n";
        }

        if($f->has_tag('gene'))
        {
          $gene = join(',',sort $f->each_tag_value('gene'));
        }
        else
        {
          if($f->has_tag('product'))
          {
            $gene = join(',',sort $f->each_tag_value('product'));
            if(length($gene)>20){ $gene = substr($gene,0,20).'..' } # avoid long product names
          }
        }

        next if($protsequence eq '');

        # NA is assigned when a valid DNA region was not reported
        $genelength= length($CDSseq) || 'NA';
        $start = $f->start(); # 1..n naturals
        $end   = $f->end();
        $strand = $f->location()->strand();
        $header = $gi." |".$taxon."|".$strain."|".basename($infile)."|".$gene."|".$genelength."|".
          "$gbaccession($source):$start-$end:$strand" .
          " ^$crossrefs^ ".$CDScoords->desc();
        if(!$already_seen{$header})
        {
          push(@{$genes{$gbaccession}},[$gi,$gene,$strand,$header,$CDSseq,$protsequence]);
          $already_seen{$header} = 1;
          $n_of_CDS++;
        }
      }
    }
  }

  # first create tmp files, in case sub does not finish cleanly
  open(DNA,">$out_dna_file.tmp") || die "# extract_CDSs_from_genbank : cannot create $out_dna_file.tmp\n";
  open(PROT,">$out_prot_file.tmp") || die "# extract_CDSs_from_genbank : cannot create $out_prot_file.tmp\n";
  foreach $gbaccession (keys(%genes))
  {
    my ($contig_n_of_genes,$contig) = (scalar(@{$genes{$gbaccession}}),$genes{$gbaccession});
    for(my $g=0;$g<scalar(@$contig);$g++)
    {
      $header = "$contig->[$g][3]|";
      if($g==0){ $header .= "neighbours:start(),"; }
      else{ $header .= "neighbours:$contig->[$g-1][0]($contig->[$g-1][2]),"; }
      if($g<$contig_n_of_genes-1)
      {
        $header .= "$contig->[$g+1][0]($contig->[$g+1][2])|";
      }
      else{ $header .= "end()|" }
      if($g==0){ $header .= "neighbour_genes:start(),"; }
      else{ $header .= "neighbour_genes:$contig->[$g-1][1],"; }
      if($g<$contig_n_of_genes-1)
      {
        $header .= "$contig->[$g+1][1]|";
      }
      else{ $header .= "end()|" }

      if($verbose)
      {
        if($contig->[$g][4] eq '')
        {
          print "# extract_CDSs_from_genbank : cannot extract DNA sequence $contig->[$g][3] from $infile\n";
        }
        elsif($contig->[$g][5] eq '')
        {
          print "# extract_CDSs_from_genbank : cannot extract protein sequence $contig->[$g][3] from $infile\n";
        }
      }

      print DNA ">$header\n$contig->[$g][4]\n";
      print PROT ">$header\n$contig->[$g][5]\n";
    }
  }
  close(DNA);
  close(PROT);

  # make final outfiles only when process finishes cleanly
  rename( "$out_prot_file.tmp" , $out_prot_file );
  rename( "$out_dna_file.tmp" , $out_dna_file );
  return $n_of_CDS;
}

# To to be used when extract_CDSs fails
sub extract_genes_from_genbank
{
  # takes a genbank input file and creates two FASTA files containing all genes in
  # dna and translated aminoacid sequences
  # returns number of genes found

  my ($infile,$out_prot_file,$out_dna_file,$verbose) = @_;
  my ($gene,$gi,$gbaccession,$header,$genelength,$CDScoords,$CDSseq,$rev);
  my ($start,$end,$strand,$source,$protsequence,$magic,$in,%already_seen,%genes);
  my $n_of_genes = 0;
  my ($taxon,$strain)=('','');
  
  if($verbose){ print "# extract_genes_from_genbank : extracting file $infile...\n" }
  
  # check input file and open it accordingly
  open(INFILE,$infile) || die "# read_FASTA_sequence_array: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    $in = new Bio::SeqIO(-file => "gzip -dc $infile |", -format => 'genbank' )
  }
  elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    $in = new Bio::SeqIO(-file => "bzip2 -dc $infile |", -format => 'genbank' )
  }
  else{ $in = new Bio::SeqIO(-file => $infile, -format => 'genbank' ) }
  
  my $seqobj = Bio::Seq->new( -display_id => 'tmp' );
  while( my $seq = $in->next_seq) # scan all sequences found in $input
  {
    $seq->alphabet('dna'); 
    $gbaccession = $seq->accession();
    if($gbaccession eq 'unknown'){ $gbaccession = $seq->display_id() } # prokka-compatible
    $source = '';
    foreach my $f ($seq->get_SeqFeatures)
    {
      if($f->primary_tag() =~ /source/)
      {
        $source = $f->end();
        if($taxon eq '' && $f->has_tag('organism'))
        {
          foreach my $element ($f->each_tag_value('organism')){ $taxon .= "[$element],"; }
          chop $taxon;
        }
        if($strain eq '' && $f->has_tag('strain'))
        {
          $strain = join(',',sort $f->each_tag_value('strain'));
        }
      }
      elsif($f->primary_tag() =~ /gene/)
      {
        next if $f->has_tag('pseudo'); # skip pseudogenes altogether

        $gene=$gi=$genelength=$protsequence=$CDSseq=$rev='';
        $CDScoords = $f->spliced_seq();

        if($f->location->isa('Bio::Location::SplitLocationI'))
        {
          for my $location ( $f->location->sub_Location() )
          {
            $CDSseq .= $seq->subseq($location->start,$location->end);
            if($location->strand == -1){ $rev = 1 }
          }

          if($rev)
          {
            $CDSseq =~ tr/acgtnACGTN/tgcanTGCAN/;
            $CDSseq = reverse($CDSseq);
          } #print "$CDSseq\n";
        }
        else
        {
          $CDSseq = $CDScoords->{'seq'};
        }

        # translate CDS sequence to protein
        $seqobj->seq($CDSseq);
        $protsequence = $seqobj->translate()->seq();chop($protsequence); 
        
	# check internal stop codons
        next if($protsequence =~ /\*/);

        if($gi eq '' && $f->has_tag('locus_tag'))
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('locus_tag')); #print "$gi\n";
        }
        if($f->has_tag('gene'))
        {
          $gene = join(',',sort $f->each_tag_value('gene'));
        }

        # NA is assigned when a valid DNA region was not reported
        $genelength= length($CDSseq) || 'NA';
        $start = $f->start(); # 1..n naturals
        $end   = $f->end();
        $strand = $f->location()->strand();
        $header = $gi." |".$taxon."|".$strain."|".basename($infile)."|".$gene."|".$genelength."|".
          "$gbaccession($source):$start-$end:$strand" .
          " ^^ ".$CDScoords->desc();
        if(!$already_seen{$header})
        {
          push(@{$genes{$gbaccession}},[$gi,$gene,$strand,$header,$CDSseq,$protsequence]);
          $already_seen{$header} = 1;
          $n_of_genes++;
        }
      }
    }
  }

  # first create tmp files, in case sub does not finish cleanly
  open(DNA,">$out_dna_file.tmp") || die "# extract_genes_from_genbank : cannot create $out_dna_file.tmp\n";
  open(PROT,">$out_prot_file.tmp") || die "# extract_genes_from_genbank : cannot create $out_prot_file.tmp\n";
  foreach $gbaccession (keys(%genes))
  {
    my ($contig_n_of_genes,$contig) = (scalar(@{$genes{$gbaccession}}),$genes{$gbaccession});
    for(my $g=0;$g<scalar(@$contig);$g++)
    {
      $header = "$contig->[$g][3]|";
      if($g==0){ $header .= "neighbours:start(),"; }
      else{ $header .= "neighbours:$contig->[$g-1][0]($contig->[$g-1][2]),"; }
      if($g<$contig_n_of_genes-1)
      {
        $header .= "$contig->[$g+1][0]($contig->[$g+1][2])|";
      }
      else{ $header .= "end()|" }
        
      if($g==0){ $header .= "neighbour_genes:start(),"; }
      else{ $header .= "neighbour_genes:$contig->[$g-1][1],"; }
      if($g<$contig_n_of_genes-1)
      {
        $header .= "$contig->[$g+1][1]|";
      }
      else{ $header .= "end()|" }
      
      print DNA ">$header\n$contig->[$g][4]\n";
      print PROT ">$header\n$contig->[$g][5]\n";
    }
  }
  close(DNA);
  close(PROT);

  # make final outfiles only when process finishes cleanly
  rename( "$out_prot_file.tmp" , $out_prot_file );
  rename( "$out_dna_file.tmp" , $out_dna_file );

  return $n_of_genes;
}

sub extract_gene_name
{
  # returns gene name from a FASTA header generated by extract_CDSs_from_genbank
  return (split(/\|/,$_[0]))[4];
}


# Input FASTA sequences are expected to be aligned
# takes 5 args: 
# 1) filename string, required
# 2) type of sequence (scalar), required 
# 3) blunt [min allowed block] (scalar), optional
# 4) array ref with sequence names of group A, optional
# 5) array ref with sequence names of group B, only if 4) is set
# 6) verbose (boolean), optional
#
# returns 7 scalars: 
# 1) ref to 2D array with 2ary indexes: NAME,SEQ; first sequence is '0'
# 2) number of SNP positions (integer)
# 3) number of pars positions (integer)
# 4) number of private (A vs B) variants (integer)
# 5) string with comma-sep SNP positions
# 6) string with comma-sep SNP positions
# 7) string with comma-sep SNP positions
# 8) number of unaligned groupA sequences (integer)
# 9) number of unaligned groupB sequences (integer)
sub check_variants_FASTA_alignment
{	
  my ( $infile, $peptideOK, $blunt, $arrayA, $arrayB, $verbose ) = @_;
  
  my (@FASTA,@matrix,%length,@groupA,@groupB);
  my ($name,$seq,$nt,$l,$seqname,$nameOK);
  my ($n_of_sequences,$missA,$missB) = (-1,0,0);

  # read in multiFASTA
  open(FASTA,"<$infile") || die "# check_variants_FASTA_alignment: cannot read $infile $!\n";
  while(<FASTA>)
  {
    next if(/^$/ || /^#/);
    if(/^\>(.*?)[\n\r]/)
    {
      $n_of_sequences++; # first sequence ID is 0
      $name = $1;
      $nt = 0;
      
      # match taxon names if arrayA/B are not null
      $nameOK = 0;
      foreach $seqname (@$arrayA)
      {
        if($name =~ m/\Q$seqname\E/)
        { 
          warn "# check_variants_FASTA_alignment: $n_of_sequences groupA\n" if($verbose);
          push(@groupA,$n_of_sequences); 
          $nameOK = 1; 
          last;
        }
      }
      
      if(!$nameOK)
      {
        foreach $seqname (@$arrayB)
        {
          if($name =~ m/\Q$seqname\E/)
          { 
            warn "# check_variants_FASTA_alignment: $n_of_sequences groupB\n" if($verbose);
            push(@groupB,$n_of_sequences); 
            $nameOK = 1; 
            last; 
          }
        }
      } 
      
      if($verbose && !$nameOK){ warn "# check_variants_FASTA_alignment: cannot match $name\n" }    
      
      # store parsed name in matrix
      $FASTA[$n_of_sequences][NAME] = $name;
    }
    elsif($n_of_sequences>-1)
    {
      $seq = $_;
      $seq =~ s/[\s|\n]//g;
      $seq = uc($seq); # does not conserve original case
      $FASTA[$n_of_sequences][SEQ] .= $seq; 
      
      # split sequence in order to later check for variants
      foreach my $letter (split(//,$seq))
      {
        $matrix[$n_of_sequences][$nt++] = $letter;
      }
    }	
  }
  close(FASTA);

  if($n_of_sequences == -1)
  {
    warn "# check_variants_FASTA_alignment: $infile contains no sequences, skip it\n";
	  return ([],-1,-1,-1,'','','',-1,-1);
  }

  # check sequences are aligned
  foreach $seq (0 .. $n_of_sequences)
  {
    $l = length($FASTA[$seq][SEQ]);
    $length{$l} = 1;
  }
	 
  if(scalar(keys(%length)) > 1)
  {
    warn "# check_variants_FASTA_alignment: input seqs at $infile are not aligned, skip it\n";
    return ([],-1,-1,-1,'','','',-1,-1);
  }
  
  # check included taxa names if required
  if(@$arrayA)
  {
    if($verbose)
    {
      printf(STDERR "# check_variants_FASTA_alignment: found %d sequences of groupA\n",scalar(@groupA));
      printf(STDERR "# check_variants_FASTA_alignment: found %d sequences of groupB\n",scalar(@groupB));
    }

    $missA = scalar(@$arrayA)-scalar(@groupA);
    $missB = scalar(@$arrayB)-scalar(@groupB);
  }
  
  # produce a block of blunt-end sequences
  my ($first_blunt,$last_blunt,$gaps) = (0,$l-1);

  if($blunt)
  {
    # first check left/5' side
    $nt = 0;
    while($nt<$last_blunt)
    {
      $gaps=0;										
      foreach $seq (0 .. $n_of_sequences) 
      {
        if($matrix[$seq][$nt] eq '-')
        {
          $gaps++;
          last;
        }
      }
      if($gaps>0)
      {
        $first_blunt++;
        $nt++;
        next;
      }
      else{ last }
    }

    # now check right/3' end
    $nt = $last_blunt;
    while($nt>$first_blunt)
    {   
      $gaps=0;                           
      foreach $seq (0 .. $n_of_sequences)
      {
        if($matrix[$seq][$nt] eq '-')
        {
          $gaps++;
          last;
        }
      }
      if($gaps>0)
      {
        $last_blunt--;
        $nt--;
        next;
      }
      else{ last }
    }

    # check length of trimmed block
    $l = ($last_blunt-$first_blunt)+1;
    if($l < $blunt)
    {
      warn "# check_variants_FASTA_alignment: trimmed block at $infile is too short ($l<$blunt), skip it\n";
      return ([],-1,-1,-1,'','','',-1,-1);
    }
  } #warn "# $first_blunt $last_blunt\n";
  
  # annotate variants within possibly trimmed block of aligned sequences 
  my @aminoacids  = qw( C S T P A G N D E Q H R K M I L V F Y W X - );
  my @nucleotides = qw( A C G T N - );
  my ($letter,$degen,$pars,$G,$snps,$letterA,$letterB);
  my ($snps_string,$pars_string,$private_string) = ('','','');
  my ($snps_total,$pars_total,$private_total) = (0,0,0);
  my (%freq,%freqA,%freqB,@trimmed_seq,@alphabet);
  my $effective_pos = 1;
  
  if($peptideOK)
  { 
    @alphabet = @aminoacids; 
    $degen = 'X';
  }
  else
  { 
    @alphabet = @nucleotides; 
    $degen = 'N';
  }

  $nt = $first_blunt;
  while($nt<=$last_blunt)
  {
    # init
    $pars=$G=$snps=0;
    foreach $letter (@alphabet){ $freq{$letter} = 0 }
    
    foreach $seq (0 .. $n_of_sequences)
    {
      $letter = $matrix[$seq][$nt]; 
      $freq{$letter}++; 
      $trimmed_seq[$seq] .= $matrix[$seq][$nt];
    }
    
    # check alignment column
    foreach $letter (@alphabet)
    {
      next if($freq{$letter} == 0 || $letter eq $degen);
      if($freq{$letter}){ $snps++ }
      if($freq{$letter} > 1){ $pars++ }
      if($letter eq '-'){ $G = 1 } # columns with gaps
    }
    
    # check SNPs
    if($snps > 1 && $G == 0)
    { 
      $snps_string .= "$effective_pos,";
      $snps_total++;
    }
    
    # parsimony-informative positions: 2+ letters observed twice
    if($pars > 1)
    {
      $pars_string .= "$effective_pos,";
      $pars_total++; 
    }
    
    # check private variants in group A vs groupB
    if(@groupA && ($snps || $G))
    {
      %freqA = ();
      %freqB = ();
      
      foreach $seq (@groupA)
      {
        $letterA = $matrix[$seq][$nt];
        $freqA{$letterA}++; 
      }
      
      foreach $seq (@groupB)
      {
        $letterB = $matrix[$seq][$nt];
        $freqB{$letterB}++; 
      } 
        
      if(scalar(keys(%freqA))== 1 && scalar(keys(%freqB))== 1 && 
        $letterA ne $degen && $letterB ne $degen && $letterA ne $letterB)
      {
        $private_string .= "$effective_pos($letterA:$letterB),";
        $private_total++;
      }
    }  
    
    $effective_pos++;
    $nt++;
  }
  
  # save possibly trimmed sequences
  foreach $seq (0 .. $n_of_sequences)
  {
    $FASTA[$seq][SEQ] = $trimmed_seq[$seq];		
  }

  return (\@FASTA,$snps_total,$pars_total,$private_total,$snps_string,$pars_string,$private_string,$missA,$missB);
}    

# Input FASTA sequences are expected to be aligned
# takes 5 args: filename (string),  type of sequence (scalar), 
# overlap (natural), maxmismatches (natural), maxgaps (natural)
# returns: array of aligned, collapsed sequences; 
# 2D array has 2ary indexes: NAME,SEQ; first sequence is '0'
sub collapse_taxon_alignments
{
  my ( $infile, $peptideOK, $min_overlap, $max_mismatches, $max_gaps, $verbose ) = @_;
  
  sub _get_mismatches_gaps
  {
    my ($ref_MSA, $seqA, $seqB, $start, $end) = @_;
    my ($mismatches,$gaps,$col,$A,$B) = (0,0);
    for $col ($start .. $end)
    {
      $A = $ref_MSA->[$seqA][$col];
      $B = $ref_MSA->[$seqB][$col];
      if($A eq $B){ next }
      elsif($A eq '-' || $B eq '-'){ $gaps++ }
      else{ $mismatches++ }
    }
    return ($mismatches,$gaps);
  }
  
  my (@FASTA,@FASTAcons,@matrix,%length,%taxon2id);
  my (%left,%right,%midpoint);
	my ($name,$seq,$seq2,$nt,$l,$taxon,$MSAwidth);
  my ($overlap,$mismatches,$gaps,$cons_start,$cons_end);
	my ($n_of_sequences,$n_of_collapsed) = (-1,0);

	# read in multiFASTA
  open(FASTA,"<",$infile) || die "# collapse_taxon_alignments: cannot read $infile $!\n";
  while(<FASTA>)
  {
		next if(/^$/ || /^#/);
    if(/^>(\S+)/)
    {
			$n_of_sequences++; # first sequence ID is 0
      $name = $1;
			$nt = 0;
      $FASTA[$n_of_sequences][NAME] = $name;
      
      $taxon = '';
      if($name =~ m/^(\S+).*?\[(\S+?)\]*$/ || $name =~ m/^(\S+).*?\[(\S+?)\]\s\|/ ||
        $name =~ m/^(\S+).*?\[(.+?)\]*\s\|/ || $name =~ m/^(\S+).*?\[(.+?)\]*/)
      {
        $taxon = (split(/\.f/,$2))[0];
        push(@{$taxon2id{$taxon}},$n_of_sequences); #print "$taxon $n_of_sequences\n";
      }  
    }
    elsif($n_of_sequences>-1)
    {
		  $seq = $_; 
		  $seq =~ s/[\s|\n]//g;
		  $seq = uc($seq); # does not conserve original case
      $FASTA[$n_of_sequences][SEQ] .= $seq; 
      
      # split sequence in order to later check for variants
	    foreach my $letter (split(//,$seq))
		  {
		    $matrix[$n_of_sequences][$nt++] = $letter;
      }  
    }  
	}
  close(FASTA);
  
	if($n_of_sequences == -1)
	{
    warn "# collapse_taxon_alignments: $infile contains no sequences, exit\n";
	  return [];
  }

  if(scalar(keys(%taxon2id)) < 1)
	{
    warn "# collapse_taxon_alignments: cannot parse taxon names in $infile, exit\n";
	  return [];
  }
  
	# check sequences are aligned
  $MSAwidth = 0;
	foreach $seq (0 .. $n_of_sequences)
	{
		$l = length($FASTA[$seq][SEQ]);
		$length{$l} = 1;
    $MSAwidth = $l;
	}
	 
	if(scalar(keys(%length)) > 1)
	{
		warn "# collapse_taxon_alignments: input seqs at $infile are not aligned, exit\n";
		return [];
	}#else{ print "# MSAwidth=$MSAwidth\n" }
  
  # compute coordinates of aligned fragments
  foreach $seq (0 .. $n_of_sequences)
  {
    # N/5' end
    $nt=0;
    while($matrix[$seq][$nt] eq '-'){ $nt++ }
    $left{$seq} = $nt;
    
    # C/3' end a d midpoint
    $nt = $MSAwidth-1;
    while($matrix[$seq][$nt] eq '-'){ $nt-- }
    $right{$seq} = $nt;
    
    $midpoint{$seq} = sprintf("%1.0f",$left{$seq} + (($right{$seq}-$left{$seq}+1)/2));  
  }    
  
  # compare fragments of same taxon
  foreach $taxon (sort(keys(%taxon2id)))
  {
    my ($merged,%collapsed);
   
    # ensure that                left-most            and   longest fragments are handled first
    my @sorted_fragments = sort {$left{$a}<=>$left{$b} || $midpoint{$b}<=>$midpoint{$a}} @{$taxon2id{$taxon}};

    warn "# collapse_taxon_alignments: $taxon (n=".scalar(@sorted_fragments).")\n"; #,".join(',',@sorted_fragments)."\n";
    
    foreach $seq (0 .. $#sorted_fragments)
    {
      next if($collapsed{$sorted_fragments[$seq]});
    
      my @consensus = @{$matrix[$sorted_fragments[$seq]]};
      
      $merged = '';
      $cons_start = $left{$sorted_fragments[$seq]};
      $cons_end   = $right{$sorted_fragments[$seq]};
      print "# $taxon $sorted_fragments[$seq] $cons_start $cons_end\n" if($verbose);
           
      foreach $seq2 ($seq+1 .. $#sorted_fragments)
      {
        next if($collapsed{$sorted_fragments[$seq2]});
        print "## $taxon $sorted_fragments[$seq2] $left{$sorted_fragments[$seq2]} $right{$sorted_fragments[$seq2]}\n" if($verbose);
        
        # $seq includes $seq2, should not happen the other way aorund as they are sorted above
        if($left{$sorted_fragments[$seq]} <= $left{$sorted_fragments[$seq2]} 
          && $right{$sorted_fragments[$seq]} >= $right{$sorted_fragments[$seq2]})
        {
          ($mismatches,$gaps) = _get_mismatches_gaps( \@matrix, $sorted_fragments[$seq], $sorted_fragments[$seq2], 
            $left{$sorted_fragments[$seq2]}, $right{$sorted_fragments[$seq2]} ); #print "seq{seq2} $mismatches $gaps\n"; 
            
          if($mismatches <= $max_mismatches && $gaps <= $max_gaps)
          {
            $collapsed{$sorted_fragments[$seq2]} = 1; 
            $merged .= "{$sorted_fragments[$seq2]},"; 
          }
        }
        #elsif($left{$sorted_fragments[$seq2]} <= $left{$sorted_fragments[$seq]} # $seq2 includes $seq
        #  && $right{$sorted_fragments[$seq2]} >= $right{$sorted_fragments[$seq]}){
        #  #  ($mismatches,$gaps) = _get_mismatches_gaps( \@matrix, $sorted_fragments[$seq], $sorted_fragments[$seq2], 
        #  #    $left{$sorted_fragments[$seq]}, $right{$sorted_fragments[$seq]} ); #print "seq2[seq] $mismatches $gaps\n"; 
        #}
        else # check end-to-start overlaps
        {        
          $overlap = ($right{$sorted_fragments[$seq]}-$left{$sorted_fragments[$seq2]})+1; 
          if($overlap > $min_overlap)
          {
            ($mismatches,$gaps) = _get_mismatches_gaps( \@matrix, $sorted_fragments[$seq], $sorted_fragments[$seq2], 
            $left{$sorted_fragments[$seq2]}, $right{$sorted_fragments[$seq]} ); #print "seq>seq2 $overlap $mismatches $gaps\n";
            
            if($mismatches <= $max_mismatches && $gaps <= $max_gaps)
            {
              $collapsed{$sorted_fragments[$seq2]} = 1;
              $merged  .= ">$sorted_fragments[$seq2],"; 
              $cons_end = $right{$sorted_fragments[$seq2]};
              
              foreach $nt ($right{$sorted_fragments[$seq]}+1 .. $cons_end)
              {
                $consensus[$nt] = $matrix[$sorted_fragments[$seq2]][$nt];
              }
            }
          }
        }
      }
      
      # add this raw/merged sequence
      $FASTAcons[$n_of_collapsed][NAME]  = $FASTA[$sorted_fragments[$seq]][NAME];
      if($merged ne '')
      {
        $FASTAcons[$n_of_collapsed][NAME] .= " collapsed:$sorted_fragments[$seq],$merged";
      }
      $FASTAcons[$n_of_collapsed][SEQ] = join('',@consensus);
      $n_of_collapsed++;
    }  
              
  } # foreach taxon
  
  return \@FASTAcons;
}

# based on http://www.unix.org.ua/orelly/perl/cookbook/ch04_18.htm
# generates a random permutation of @array in place and returns 
# string with concatenated elements
sub fisher_yates_shuffle
{
  my ($array) = (@_);
  my ($i,$j);

  for ($i = @$array; --$i; )
  {
    $j = int(rand($i+1));
    next if $i == $j;
    @$array[$i,$j] = @$array[$j,$i];
  }

  return join('',@$array);
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

# works only on Linux/macOS systems, which have 'ps'
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

# called from get_homologues*.pl to terminate warning that a software requirement is not met
# uses globals: $0
sub warn_missing_soft
{
  my ($soft) = @_;
  die "# ERROR: cannot run $soft as required software is not in place\n".
  	"# EXIT : run $0 -v to check all required binares and data sources\n";
}



1;
