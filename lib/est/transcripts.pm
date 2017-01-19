package transcripts;

# Code library created by Bruno Contreras-Moreira (2015-17)
# mainly for transcripts2cds.pl

use strict;
require Exporter;
use File::Temp;
use phyTools;
use File::Basename;

our @ISA = qw( Exporter );

our @EXPORT = qw(

  shorten_headers_FASTA_file
  executeMAKEDB
  format_DIAMOND_command
  executeFORMATDB_EST
  format_BLASTX_command
  parse_blastx_cds_sequences
  parse_transcoder_sequences
  sequence_consensus
  is_compressed
  extract_compressed_file
  executeGMAPBUILD
  format_GMAP_command
  parse_GMAP_results
  format_transdecoder_command

  $BLASTXDB
  $BLASTX_NOCPU
  $MINORFLENGTH
  $MAXBLASTXEVALUE
  $MINCONOVERLAP

  );

# all these defined in phyTools.pm:
our $BLASTX_NOCPU        = 2; # most current cpus are multicore
our $MINORFLENGTH        = 50;
our $MAXBLASTXEVALUE     = 1e-5; # as in http://proteomics.ysu.edu/tools/OrfPredictor.html
our $BLASTX              = $ENV{"EXE_BLASTX_EST"};
our $FORMATDB            = $ENV{"EXE_FORMATDB_EST"};
our $TRANSDECODEXE       = $ENV{"EXE_TRANSDECOD_EST"};
our $MAKEDB              = $ENV{"EXE_DMNFT"};
our $DIAMONDEXE          = $ENV{"EXE_DMNDX_EST"};  
our $BLASTXDB            = $ENV{"BLASTXDB"};
our $GMAPEXE             = $ENV{"EXE_GMAP"};
our $GMAPBUILDEXE        = $ENV{"EXE_GMAPBUILD"};

# runtime & cutoff variables
my $MAXBLASTXHITS        = 5;
my $MINCONOVERLAP        = 90; # number of bases to call an overlap while calculating a consensus (30 aminos)
my $MAXACCESSIONLENGTH   = 30;

# gets a possibly compressed FASTA file and returns the name of a new one with short FASTA headers,
# Transdecoder does not like them, and upper case sequences
# Updated Jan2017
sub shorten_headers_FASTA_file
{
  my ($infile,$outfile) = @_;

  my $fasta_ref = read_FASTA_file_array($infile);
  
  open(SHORT,'>',$outfile) || 
    die "# short_headers_FASTA_file: cannot create $outfile\n";
  
  my $n_of_seqs = 0;
  foreach my $seq ( 0 .. $#{$fasta_ref} )
  {
    my ($acc,$rest) = split(/\s+/,$fasta_ref->[$seq][NAME],2);
    if(length($acc) > $MAXACCESSIONLENGTH)
    {
      ($acc,$rest) = split(/[\|\[]/,$fasta_ref->[$seq][NAME],2);
    }
    
    printf( SHORT ">%s %s\n%s\n",$acc,$rest,uc($fasta_ref->[$seq][SEQ]));
    $n_of_seqs++;
  }
  
  close(SHORT);
  
  unlink($outfile) if(!$n_of_seqs);
  
  return $n_of_seqs;
}


sub is_compressed
{
  my ($infile) = @_;
  my ($magic);

  open(INFILE,$infile) || die "# is_compressed: cannot open $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    return 'gzip';
  }
  elsif($infile =~ /\.bz2/ || $magic eq "BZ") # BZIP2 compressed input
  {
    return 'bzip2';
  }

  return 0;
}

# returns name of temporary file with extracted content of $infile
sub extract_compressed_file
{
  my ($infile,$format) = @_;

  # create temporary filename
  my $tmpfilename = tmpnam();

  if($format eq 'gzip')
  {
    system("gzip -dc $infile > $tmpfilename");
    return $tmpfilename;
  }
  elsif($format eq 'bzip2')
  {
    system("bzip2 -dc $infile > $tmpfilename");
    return $tmpfilename;
  }
  else
  {
    warn "# extract_compressed_file: unrecognized format ($format), exit\n";
    return '';
  }
}

sub executeMAKEDB
{
  my ($in) = @_;

  my ($command);

  print("\n# running makedb with $in\n");
  if(-s $in)
  {
    $command = "$MAKEDB --in $in --db $in"; 

    open(EXE,"$command |") || die "# ERROR (executeMAKEDB) : cannot run $command : $!\n";
    while(<EXE>){}
    close(EXE);
  }
  else{ die "# executeMAKEDB : cannot find input FASTA file $in\n"; }
}

sub format_DIAMOND_command
{
  my ($infile,$outfile,$db,$Evalue,$gencode,$plustStrandOnly,$hits_to_show) = @_;
  
  my $command = "$DIAMONDEXE -p $BLASTX_NOCPU -q $infile --evalue $Evalue -d $db -o $outfile " .
    "--max-target-seqs $MAXBLASTXHITS --quiet --more-sensitive ";

  if(defined($gencode) && $gencode > 1){ $command .= " --query-gencode $gencode " }
  if($plustStrandOnly){ $command .= ' --forwardonly ' }
  
  $command .= '--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ';

  return $command;
}

sub executeFORMATDB_EST
{
  my ($in) = @_;

  my ($command);

  print("\n# running makeblastdb with $in\n");
  if(-s $in)
  {
    $command = "$FORMATDB -in $in -dbtype prot ";

    open(EXE,"$command |") || die "# ERROR (executeFORMATDB_EST) : cannot run $command : $!\n";
    while(<EXE>){}
    close(EXE);
  }
  else{ die "# executeFORMATDB_EST : cannot find input FASTA file $in\n"; }
}

sub format_BLASTX_command
{
  my ($infile,$outfile,$db,$Evalue,$gencode,$plustStrandOnly,$hits_to_show) = @_;

  #my $command = "$BLASTX -num_threads $BLASTX_NOCPU -outfmt \"6 std qlen qseq sseq\" " . #debug
  my $command = "$BLASTX -num_threads $BLASTX_NOCPU -outfmt \"6 std qlen\" " .
    "-query $infile -evalue $Evalue -db $db -out $outfile " .
    "-max_target_seqs $MAXBLASTXHITS ";

  if(defined($gencode) && $gencode > 1){ $command .= " -query_gencode $gencode " }
  if($plustStrandOnly){ $command .= ' -strand plus ' }

  return $command;
}

# parses longest ORF detected by transdecoder
# Updated Jan2017
sub parse_transcoder_sequences
{
  my ($infile) = @_;

  my $fasta_ref = read_FASTA_file_array($infile);
  my ($seqname,%seqs);

  foreach my $seq ( 0 .. $#{$fasta_ref} )
  {
    #>AT1G02065.2|m.737 AT1G02065.2|g.737 type:complete len:247 AT1G02065.2:59-799(+)
    #>AT1G02065.2|m.738 AT1G02065.2|g.738 type:5prime_partial len:168 AT1G02065.2:3-506(+)
    #>AT1G02065.2|m.739 AT1G02065.2|g.739 type:complete len:78 AT1G02065.2:1474-1241(-)
    #>AT1G02065.2|m.740 AT1G02065.2|g.740 type:complete len:78 AT1G02065.2:373-140(-)
    #>AT1G02065.2|m.741 AT1G02065.2|g.741 type:complete len:75 AT1G02065.2:707-483(-)
    #>AT1G02065.2|m.742 AT1G02065.2|g.742 type:complete len:50 AT1G02065.2:1220-1071(-)

    $seqname = (split(/\|[mg]\./,$fasta_ref->[$seq][NAME]))[0];

    if(!$seqs{$seqname} || # parse only longest ORF
      length($fasta_ref->[$seq][SEQ]) > length($seqs{$seqname}))
    {
      $seqs{$seqname} = $fasta_ref->[$seq][SEQ];
      #print "$seqname $fasta_ref->[$seq][SEQ]\n";
    }
  }

  return \%seqs;
}

# parses top matches (on the same frame) of each query
# updated Jan2017
sub parse_blastx_cds_sequences
{
  my ($FASTAdnafile,$blastxfile) = @_;

  my ($qseqid,$sseqid,$qstart,$qend,$qlen,$qseq,$frame);
  my ($magic,$offset5,$offset3,$pos,$nt,$revcomp,$lastframe,$lastrevcomp);
  my (%blastxseqs,%blastxhits,%blastxrev,%seqs,@alnseq);

  my $fasta_ref = read_FASTA_file_array($FASTAdnafile);
  my ($seqname,%dnaseqs);

  foreach my $seq ( 0 .. $#{$fasta_ref} )
  {
    $seqname = (split(/\s+/,$fasta_ref->[$seq][NAME]))[0];
    $dnaseqs{$seqname} = uc($fasta_ref->[$seq][SEQ]);
    #print "$seqname $fasta_ref->[$seq][SEQ]\n";
  }

  # check input file format and open it accordingly
  open(INFILE,$blastxfile) || die "# parse_blastx_cds_sequences: cannot read $blastxfile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($blastxfile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    if(!open(BLASTX,"gzip -dc $blastxfile |"))
    {
      die "# parse_blastx_cds_sequences: cannot read GZIP compressed $blastxfile $!\n"
        ."# please check gzip in installed\n";
    }
  }
  elsif($blastxfile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    if(!open(BLASTX,"bzip2 -dc $blastxfile |"))
    {
      die "# parse_blastx_cds_sequences: cannot read BZIP2 compressed $blastxfile $!\n"
        ."# please check bzip2 in installed\n";
    }
  }
  else{ open(BLASTX,$blastxfile) || die "# parse_blastx_cds_sequences: cannot read $blastxfile $!\n"; }

  while(<BLASTX>)
  {
    #($qseqid,$sseqid,$pident,$length,$mismatch,$gaps,$qstart,$qend,$sstart,$send,$evalue,$score,$qlen) = split('\t',$_);
    #comp37170_c0_seq1   gi|514813219|ref|XP_004981393.1|    41.49   188 76  6   941 384 1   156 3e-31   ' 116'    960
    #s/ //g;# bitscore is space-padded at least on ncbi-blast-2.2.27+
    #if(/^(\S+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t(\S+)/)
    if(/^(\S+)\t(\S+)\t\S+\t\S+\t\S+\t\S+\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t\s*\S+\t(\S+)/)
    {
      ($qseqid,$sseqid,$qstart,$qend,$qlen) = ($1,$2,$3,$4,$5);
      #print "$qseqid,$sseqid,$qstart,$qend,$qlen\n";

      if($qstart>$qend)
      {
        $revcomp = 1;
        ($qstart,$qend) = ($qend,$qstart);
      }
      else{ $revcomp = 0 }

      $frame = $qstart % 3;
      if($revcomp){ $frame = -$frame }

      # 1st hit for this query, initialize sequence array with lower case nucleotides
      if(!$blastxseqs{$qseqid})
      {
        @alnseq = split(//,lc($dnaseqs{$qseqid}));
        $pos = 0;
        foreach $nt (@alnseq){ $blastxseqs{$qseqid}[$pos++] = $nt };
        $blastxhits{$qseqid} = $sseqid;
        $blastxrev{$qseqid} = $revcomp;
        $lastframe = $frame;
        $lastrevcomp = $revcomp;
      }
      else
      {
        if($revcomp != $lastrevcomp)
        {
          warn "# parse_blastx_cds_sequences: likely chimera: $qseqid\n";
        }
        
        # take only additional HSPs in the same frame
        next if($frame != $lastframe);
      } 

      # add aligned sequence stretch (upper case)
      $qseq = substr($dnaseqs{$qseqid},$qstart-1,$qend-$qstart+1); 
      @alnseq = split(//,$qseq);
      $pos = $qstart-1;
      foreach $nt (@alnseq){ $blastxseqs{$qseqid}[$pos++] = $nt }
    }
  }
  close(BLASTX);

  foreach $qseqid (keys(%blastxseqs))
  {
    $qseq = join('',@{$blastxseqs{$qseqid}});
    
    # remove non-matched termini
    $qseq =~ s/^[a-z]+//g;
    $qseq =~ s/[a-z]+$//g; 
    
    # get reverse complement of negative strand alignments
    if($blastxrev{$qseqid}){ $qseq =~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/; $qseq = reverse($qseq) }
    
    $seqs{$qseqid} = uc($qseq); #print ">$qseqid\n$qseq\n";
  }

  return (\%seqs,\%blastxhits);
}

# Expects upper case nucleotide sequences
# returns 1) consensus sequence and 2) attached evidence
# Updated Jan2017
sub sequence_consensus
{
  my ($seq1,$seq2,$label1,$label2,$priority,$verbose) = @_;

  my $ref_alignments = local_alignments($seq1,$seq2);
  
  my @align = @{$ref_alignments->[0]}; # take only first alignment

  #print "$align[0] ".
  #	"$align[1] $align[2] $align[3] ". # 1
  #	"$align[4] $align[5] $align[6] $align[7]\n"; # 2

  if($align[0] >= $MINCONOVERLAP)
  {
    if($align[2] == $align[3] && $align[4] == 1)
    {
      # 1----------
      #      2-----------
      #print "# sequence_consensus: 1 ($align[0])\n" if($verbose);
      #print "$seq1\n$seq2\n$align[7]\n";      
      return ( substr($seq1,0,$align[1]-1).$seq2 , "$label1.$label2");
    }
    elsif($align[1] == 1 && $align[5] == $align[6])
    {
      #     1----------
      # 2-----------
      #print "# sequence_consensus: 2 ($align[0])\n" if($verbose);
      #print "$seq1\n$seq2\n$align[7]\n";      
      return ( $seq2.substr($seq1,$align[2]) , "$label2.$label1");
    }
    elsif($align[4] == 1 && $align[5] == $align[6])
    {
      # 1--------------------
      #      2-----------
      #print "# sequence_consensus: 3 ($align[0])\n" if($verbose);
      return ( $seq1 , "$label1<$label2" );
    }
    elsif($align[1] == 1 && $align[2] == $align[3])
    {
      #     1----------
      # 2------------------
      #print "# sequence_consensus: 4 ($align[0])\n" if($verbose);
      return ( $seq2 , "$label2<$label1" );
    }
    else
    {
      if($priority == 1)
      {
        #print "# sequence_consensus: 5 ($align[0])\n" if($verbose);
        return ( $seq1 , "$label1-mismatches" );
      }
      else
      {
        #print "# sequence_consensus: 6 ($align[0])\n" if($verbose);
        #print "$seq1\n$seq2\n$align[7]\n";      
        return ( $seq2 , "$label2-mismatches" );
      }
    }
  }
  else
  {
    if($priority == 1)
    {
      return ( $seq1 , "$label1-noover" );

      #print "# sequence_consensus: 7 ($align[0])\n" if($verbose);
    }
    else
    {
      return ( $seq2 , "$label2-noover" );

      #print "# sequence_consensus: 8 ($align[0])\n" if($verbose);
    }
  }
}

sub local_alignments
{
  # adapted from http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Longest_common_substring
  # assumes only $str2 can contain undetermined chars 'X'
  
  my ($str1, $str2) = @_;
  my $alignlength = 0; # length of longest common substring
  my $len1 = length $str1;
  my $len2 = length $str2;
  my @char1 = (undef, split(//, $str1)); # $str1 as array of chars, indexed from 1
  my @char2 = (undef, split(//, $str2)); # $str2 as array of chars, indexed from 1
  my @lc_suffix; # "longest common suffix" table
  my @substrings; # list of common substrings of length $l_length
  my ($n1,$n2);

  for $n1 ( 1 .. $len1 )
  {
    for $n2 ( 1 .. $len2 )
    {
      if ($char1[$n1] eq $char2[$n2] || $char1[$n1] eq "*" || $char2[$n2] eq "N")
      {
        # We have found a matching character. Is this the first matching character, or a
        # continuation of previous matching characters? If the former, then the length of
        # the previous matching portion is undefined; set to zero.
        $lc_suffix[$n1-1][$n2-1] ||= 0;

        $lc_suffix[$n1][$n2] = $lc_suffix[$n1-1][$n2-1] + 1;

        # If the resulting substring is longer than previous ...
        if ($lc_suffix[$n1][$n2] > $alignlength)
        {

          # ... we record its length as our new max length ...
          $alignlength = $lc_suffix[$n1][$n2];

          # ... and clear our result list of shorter substrings.
          @substrings = ();
        }

        # If this substring is equal to our longest ...
        if ($lc_suffix[$n1][$n2] == $alignlength)
        {

          # ... add it to our list of solutions.
          push(@substrings, [$alignlength, ($n1-$alignlength+1), $n1, $len1, ($n2-$alignlength+1),
              $n2, $len2, substr($str1,($n1-$alignlength), $alignlength) ] );
        }
      }
    }
  }

  return \@substrings;
}

sub executeGMAPBUILD
{
  my ($path,$dbname,$infile) = @_;
  my ($command);

  print("\n# running gmapbuild with $infile\n");
  if(-s $infile)
  {
    $command = "$GMAPBUILDEXE -d $dbname -D $path $infile"; #die "$command\n";

    open(EXE,"$command 2>&1 > /dev/null |") || die "# ERROR (executeGMAPBUILD) : cannot run $command : $!\n";
    while(<EXE>){}
    close(EXE);
  }
  else{ die "# executeGMAPBUILD : cannot find input FASTA file $infile\n"; }
}

sub format_GMAP_command
{
  my ($path,$dbname,$infile,$outfile,$n_of_cpus,$hits_to_show) = @_;

  my $cpus = 1;
  if($n_of_cpus && $n_of_cpus > 1){ $cpus = $n_of_cpus }

  my $command = "$GMAPEXE -S -D $path -d $dbname --no-chimeras -t $cpus "; # standard settings
  #my $command = "$GMAPEXE -A -D $path -d $dbname --no-chimeras -t $cpus "; # debug

  if($hits_to_show){ $command .= "-n $hits_to_show " }

  $command .= " $infile > $outfile  2> /dev/null";

  return $command;
}

# reads FASTA sequences and gmap report and returns:
# 1) best match coordinates, if any
# 2) cds sequences made of concatenated exons, else original sequences
sub parse_GMAP_results
{
  my ($FASTAdnafile,$gmapfile,$verbose) = @_;

  my ($magic,$seqname,$exon,%dnaseqs,%gmaphit,%cds,@intronL);
  my ($qseqid,$sseqid,$sstart,$send,$strand,$exstart,$exend);

  my $fasta_ref = read_FASTA_file_array($FASTAdnafile);
  foreach my $seq ( 0 .. $#{$fasta_ref} )
  {
    $seqname = (split(/\s+/,$fasta_ref->[$seq][NAME]))[0];
    $dnaseqs{$seqname} = $fasta_ref->[$seq][SEQ];
  }

  # check input file format and open it accordingly
  open(INFILE,$gmapfile) || die "# parse_GMAP_results: cannot read $gmapfile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($gmapfile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    if(!open(GMAP,"gzip -dc $gmapfile |"))
    {
      die "# parse_GMAP_results: cannot read GZIP compressed $gmapfile $!\n"
        ."# please check gzip in installed\n";
    }
  }
  elsif($gmapfile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    if(!open(GMAP,"bzip2 -dc $gmapfile |"))
    {
      die "# parse_GMAP_results: cannot read BZIP2 compressed $gmapfile $!\n"
        ."# please check bzip2 in installed\n";
    }
  }
  else{ open(GMAP,$gmapfile) || die "# parse_GMAP_results: cannot read $gmapfile $!\n"; }

  # add sequences matched in gmap report
  while(<GMAP>)
  {

#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore <- blast
#>952 <- gmap
#Paths (1):
#  Path 1: query 1..1719 (1719 bp) => genome 104705:546..2,167 (1622 bp)
# 	cDNA direction: indeterminate
# 	Coverage: 100.0 (query length: 1719 bp)
# 	Trimmed coverage: 100.0 (trimmed length: 1719 bp, trimmed region: 1..1719)
# 	Percent identity: 94.2 (1619 matches, 3 mismatches, 97 indels, 0 unknowns)
#Alignments:
#  Alignment for path 1:
#
#    +MLOC_67780.1:1-130  (228-357)   99% ->   ...101...  0.987, 0.989
#    +MLOC_67780.1:232-312  (358-438)   100% ->   ...1720...  0.996, 0.997
#    +MLOC_67780.1:2033-2208  (439-614)   100% ->   ...137...  0.983, 0.993
#    +MLOC_67780.1:2346-2655  (615-924)   100%
#
#  Alignment for path 2:
#
#    -MLOC_67781.1:6412-6283  (228-357)   99% ->   ...101...  0.987, 0.989
#    -MLOC_67781.1:6181-6101  (358-438)   100% ->   ...1720...  0.996, 0.997
#    -MLOC_67781.1:4380-4205  (439-614)   100% ->   ...137...  0.983, 0.993
#    -MLOC_67781.1:4067-3758  (615-924)   100%

    if(/^>(\S+)/){ $qseqid = $1 }
    elsif(/^  Path \d+: query \S+ \(\S+ bp\) => genome (\S+?):(\S+?)\.\.(\S+)/)
    {
      ($sseqid,$sstart,$send) = ($1,$2,$3);
      $sstart =~ s/,//g;
      $send =~ s/,//g;

      if($sstart>$send)
      {
        $strand = '-';
        ($sstart,$send) = ($send,$sstart);
      }
      else{ $strand = '+' }

      $gmaphit{$qseqid} = "$strand$sseqid:$sstart:$send"; #print "$strand$sseqid:$sstart:$send\n";
    }
    elsif(/^    [+-]\S+?:\d+-\d+  \((\d+)-(\d+)\)\s+\S+/)
    {
      ($exstart,$exend) = ($1,$2);

      $exon = substr($dnaseqs{$qseqid},$exstart-1,$exend-$exstart+1);
      if($strand eq '-'){ $exon =~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/; $exon = reverse($exon) }

      # calculate exon coordinates within CDS (for future use?)
      #$exstart = 1;
      #if($cds{$qseqid}){ $exstart += length($cds{$qseqid}) }
      #$exend = $exstart + length($exon) - 1;
      #$exons{$qseqid} .= "$exstart-$exend,";

      $cds{$qseqid} .= $exon;	#print "$cds{$qseqid}\n";

      # is there an intron?
      if(/   \.\.\.(\d+)\.\.\./){ push(@intronL,$1); }
    }
  }
  close(GMAP);

  # finally add sequences with no gmap matches
  foreach my $seq ( 0 .. $#{$fasta_ref} )
  {
    $seqname = (split(/\s+/,$fasta_ref->[$seq][NAME]))[0];
    next if($gmaphit{$seqname});
    $cds{$seqname} = $fasta_ref->[$seq][SEQ];
  }

  # print some stats if required
  if($verbose)
  {
    my $medianL = 0;
    my $mid = int(scalar(@intronL)/2);
    my @sortedL = sort {$a<=>$b} (@intronL);

    if(scalar(@sortedL) % 2){ $medianL = $sortedL[ $mid ] }
    else{ $medianL = ($sortedL[$mid-1] + $sortedL[$mid])/2 }

    printf("# parse_GMAP_results: matches=%d introns=%d median_intron_length=%1.0f\n",
      scalar(keys(%gmaphit)),scalar(@intronL),$medianL);
  }

  return (\%cds,\%gmaphit);
}

sub format_transdecoder_command
{
  my ($infile,$ORFlength,$plusStrand,$gencode) = @_;

  my $command = "$TRANSDECODEXE -m $ORFlength -t $infile ";
  if($plusStrand){ $command .= ' -S' }
  if($gencode > 1){ $command .= " -G $gencode" }

  return $command;
}

1;
