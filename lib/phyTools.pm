# Bruno Contreras-Moreira, Pablo Vinuesa
# 2005-16 CCG/UNAM, Mexico, EEAD/CSIC, Zaragoza, Spain
# This is a collection of subroutines used in our projects,
# including primers4clades and get_homologues.pl

package phyTools;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(
  set_phyTools_env $ERROR read_FASTA_sequence feature_is_installed check_installed_features
  extract_gene_name extract_CDSs_from_genbank extract_intergenic_from_genbank extract_features_from_genbank extract_genes_from_genbank
  find_taxa_FASTA_headers find_taxa_FASTA_array_headers read_FASTA_file_array run_PARS check_variants_FASTA_alignment NAME SEQ 
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
  'EXE_BLASTX_EST'=>'BLAST','EXE_FORMATDB'=>'database','EXE_TRANSDECOD_EST'=>'transcript',
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
  if( ! defined($ENV{"PFAMDB"}) ) { $ENV{"PFAMDB"} = $ENV{'MARFIL'}."db/Pfam-A.hmm"; } # HMMER3.0 is not compatible with Pfam > 27
  if( ! defined($ENV{'BLAST_PATH'}) ){ $ENV{'BLAST_PATH'} = $ENV{'MARFIL'}.'bin/ncbi-blast-2.2.27+/bin/'; }
  if( ! defined($ENV{'EXE_BLASTP'}) ){ $ENV{'EXE_BLASTP'} = $ENV{'BLAST_PATH'}.'blastp'; }
  if( ! defined($ENV{'EXE_BLASTN'}) ){ $ENV{'EXE_BLASTN'} = $ENV{'BLAST_PATH'}.'blastn'; }
  if( ! defined($ENV{'EXE_FORMATDB'}) ){ $ENV{'EXE_FORMATDB'} = $ENV{'BLAST_PATH'}.'makeblastdb'; }
  if( ! defined($ENV{"EXE_SPLITBLAST"}) ){ $ENV{"EXE_SPLITBLAST"} = $ENV{'MARFIL'}."_split_blast.pl"; }
  if( ! defined($ENV{"EXE_SPLITHMMPFAM"}) ){ $ENV{"EXE_SPLITHMMPFAM"} = $ENV{'MARFIL'}."_split_hmmscan.pl"; }
  if( ! defined($ENV{'EXE_PARS'}) ){ $ENV{'EXE_PARS'} = $ENV{'MARFIL'}.'bin/phylip-3.695/exe/pars'; }
  #if( ! defined($ENV{"EXE_MCL"}) ){ $ENV{"EXE_MCL"} = $ENV{'MARFIL'}."/bin/mcl-02-063/shmcl/mcl"; }
  if( ! defined($ENV{"EXE_MCL"}) ){ $ENV{"EXE_MCL"} = $ENV{'MARFIL'}."/bin/mcl-14-137/src/shmcl/mcl"; }
  if( ! defined($ENV{"EXE_INPARA"}) ){ $ENV{"EXE_INPARA"} = $ENV{'MARFIL'}."_cluster_makeInparalog.pl"; }
  if( ! defined($ENV{"EXE_ORTHO"}) ){ $ENV{"EXE_ORTHO"} = $ENV{'MARFIL'}."_cluster_makeOrtholog.pl"; }
  if( ! defined($ENV{"EXE_HOMOL"}) ){ $ENV{"EXE_HOMOL"} = $ENV{'MARFIL'}."_cluster_makeHomolog.pl"; }
  if( ! defined($ENV{"EXE_ISOFORM"}) ){ $ENV{"EXE_ISOFORM"} = $ENV{'MARFIL'}."_cluster_makeIsoform.pl"; }
  #if( ! defined($ENV{"EXE_HMMPFAM"}) ){ $ENV{"EXE_HMMPFAM"} = $ENV{'MARFIL'}."/bin/hmmer-3.0/src/hmmscan64 --noali --acc --cut_ga "; } 
  if( ! defined($ENV{"EXE_HMMPFAM"}) ){ $ENV{"EXE_HMMPFAM"} = $ENV{'MARFIL'}."/bin/hmmer-3.1b2/binaries/hmmscan --noali --acc --cut_ga "; } 
  if( ! defined($ENV{"EXE_MAKEHASH"}) ){ $ENV{"EXE_MAKEHASH"} = $ENV{'MARFIL'}."/bin/COGsoft/COGmakehash/COGmakehash "; }
  if( ! defined($ENV{"EXE_READBLAST"}) ){ $ENV{"EXE_READBLAST"} = $ENV{'MARFIL'}."/bin/COGsoft/COGreadblast/COGreadblast "; }
  if( ! defined($ENV{"EXE_COGLSE"}) ){ $ENV{"EXE_COGLSE"} = $ENV{'MARFIL'}."/bin/COGsoft/COGlse/COGlse "; }
  if( ! defined($ENV{"EXE_COGTRI"}) ){ $ENV{"EXE_COGTRI"} = $ENV{'MARFIL'}."/bin/COGsoft/COGtriangles/COGtriangles "; }
  if( ! defined($ENV{"EXE_MVIEW"}) ){ $ENV{"EXE_MVIEW"} = $ENV{'MARFIL'}."lib/mview/bin/mview "; }
  
  # transcripts/ETSs
  if( ! defined($ENV{"EXE_TRANSDECOD_EST"}) ){ $ENV{"EXE_TRANSDECOD_EST"} = $ENV{'MARFIL'}."/lib/est/TransDecoder_r20140704/TransDecoder "; }
  if( ! defined($ENV{'BLAST_PATH_EST'}) ){ $ENV{'BLAST_PATH_EST'} = $ENV{'BLAST_PATH'}; }
  if( ! defined($ENV{'EXE_BLASTX_EST'}) ){ $ENV{'EXE_BLASTX_EST'} = $ENV{'BLAST_PATH_EST'}.'blastx'; }
  if( ! defined($ENV{'EXE_FORMATDB_EST'}) ){ $ENV{'EXE_FORMATDB_EST'} = $ENV{'BLAST_PATH_EST'}.'makeblastdb'; }
  
  # diamond software
  if( ! defined($ENV{'DMND_PATH'}) ){ $ENV{'DMND_PATH'} = $ENV{'MARFIL'}.'bin/diamond-0.8.25/'; }
  if( ! defined($ENV{'EXE_DMNFT'}) ){ $ENV{'EXE_DMNFT'} = $ENV{'DMND_PATH'}.'diamond makedb'; }
  if( ! defined($ENV{'EXE_DMNDP'}) ){ $ENV{'EXE_DMNDP'} = $ENV{'DMND_PATH'}.'diamond blastp'; }
  if( ! defined($ENV{'EXE_DMNDX_EST'}) ){ $ENV{'EXE_DMNDX_EST'} = $ENV{'DMND_PATH'}.'diamond blastx'; }
  
  if( ! defined($ENV{'BLASTXDB'}) ){ $ENV{'BLASTXDB'} = $ENV{'MARFIL'}."db/uniprot_sprot.fasta"; }
  
  #if( ! defined($ENV{"EXE_GMAP"}) ){ $ENV{"EXE_GMAP"} = $ENV{'MARFIL'}."/bin/est/gmap-2014-02-20/src/gmap "; }
  #if( ! defined($ENV{"EXE_GMAPBUILD"}) ){ $ENV{"EXE_GMAPBUILD"} = $ENV{'MARFIL'}."/bin/est/gmap-2014-02-20/util/gmap_build "; }

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
      if($FASTA{$seq}{'SEQ'} !~ /[QEIPDFHKLMV]/ && length($FASTA{$seq}{'SEQ'})%3)
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
      if($FASTA[$seq][SEQ] =~ /[QEYIDFHKLVM]/i)
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

# Updated Jun2016
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
  
        if($f->has_tag('locus_tag') && $gi eq '')
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
        "neighbour_genes:$genes[$gen-1][4],$genes[$gen][4]|$taxon\n";

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
      elsif($f->has_tag('strain') && $strain eq '')
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

        if($f->has_tag('locus_tag') && $gi eq '')
        {
          $gi = "ID:".join(',',sort $f->each_tag_value('locus_tag'));
        }
        elsif($f->has_tag('db_xref') && $gi eq '')
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
            $header = $gi." |".$taxon."|".$strain."|".$genename."|".$genelength."|".
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

# Updated Oct2016
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
    $source = '';
    foreach my $f ($seq->get_SeqFeatures)
    {
      if($f->primary_tag() =~ /source/)
      {
        $source = $f->end();
        if($f->has_tag('organism') && $taxon eq '')
        {
          foreach my $element ($f->each_tag_value('organism')){ $taxon .= "[$element],"; }
          chop $taxon;
        }
        if($f->has_tag('strain') && $strain eq '')
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
          if($crossrefs =~ /SEED:fig\|([\w\.]+)/){ $gi = $1; } # Rast support Jul2016
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
        $header = $gi." |".$taxon."|".$strain."|".$gene."|".$genelength."|".
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

# Updated Oct2016
# To to be used when extract_CDSs fails, in cases such as FN869568.gbk
sub extract_genes_from_genbank
{
  # takes a genbank input file and creates two FASTA files containing all genes/CDSs in
  # aminoacid and dna sequences, respectively
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
    $source = '';
    foreach my $f ($seq->get_SeqFeatures)
    {
      if($f->primary_tag() =~ /source/)
      {
        $source = $f->end();
        if($f->has_tag('organism') && $taxon eq '')
        {
          foreach my $element ($f->each_tag_value('organism')){ $taxon .= "[$element],"; }
          chop $taxon;
        }
        if($f->has_tag('strain') && $strain eq '')
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
        $protsequence = $seqobj->translate()->seq();chop($protsequence); # faster
        #$protsequence = translate_dna2prot($CDScoords->{'seq'}); chop($protsequence);
        #print "$CDScoords->{'seq'}\n$protsequence\n";
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
        $header = $gi." |".$taxon."|".$strain."|".$gene."|".$genelength."|".
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
      if($contig->[$g-1][1] || $contig->[$g+1][1])
      {
        if($g==0){ $header .= "neighbour_genes:start(),"; }
        else{ $header .= "neighbour_genes:$contig->[$g-1][1],"; }
        if($g<$contig_n_of_genes-1)
        {
          $header .= "$contig->[$g+1][1]|";
        }
        else{ $header .= "end()|" }
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
  return $n_of_genes;
}

sub extract_gene_name
{
  # returns gene name from a FASTA header generated by extract_CDSs_from_genbank
  return (split(/\|/,$_[0]))[3];
}


# Input FASTA sequences are expected to be aligned
# takes 3 args: filename string, type of sequence (scalar) and blunt [min allowed block] (scalar)
# returns: ref to 2D array with 2ary indexes: NAME,SEQ; first sequence is '0'
sub check_variants_FASTA_alignment
{	
  my ( $infile, $peptideOK, $blunt ) = @_;
  
  my (@FASTA,@matrix,%length);
	my ($name,$seq,$nt,$l);
	my $n_of_sequences = -1;

	# read in multiFASTA
  open(FASTA,"<$infile") || die "# read_trim_FASTA_file: cannot read $infile $!\n";
  while(<FASTA>)
  {
		next if(/^$/ || /^#/);
    if(/^\>(.*?)[\n\r]/)
    {
			$n_of_sequences++; # first sequence ID is 0
      $name = $1;
			$nt = 0;
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
	  return [];
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
		return [];
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
		  return [];
	  }
	} #warn "# $first_blunt $last_blunt\n";
  
  # annotate variants within possibly trimmed block of aligned sequences 
	my @aminoacids  = qw( C S T P A G N D E Q H R K M I L V F Y W X - );
  my @nucleotides = qw( A C G T N - );
  my ($letter,$degen,$pars,$G,$snps,$snps_string,$pars_string);
  my (%freq,@trimmed_seq,@alphabet);
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
    }
    
		# parsimony-informative positions: 2+ letters observed twice
		if($pars > 1)
    {
      $pars_string .= "$effective_pos,";
    }
    
    $effective_pos++;
    $nt++;
  }
  
  # save possibly trimmed sequences
  foreach $seq (0 .. $n_of_sequences)
	{
    $FASTA[$seq][SEQ] = $trimmed_seq[$seq];		
    $FASTA[$seq][NAME] .= " SNPs:$snps_string" if($snps_string && !$peptideOK);		
  }

	return \@FASTA;
}    

1;
