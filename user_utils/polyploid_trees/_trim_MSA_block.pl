#!/usr/bin/perl -w

# This script takes a FASTA multiple alignment of [nucleotide] sequences of
# diploid and polyploid species and produces a trimmed MSA suitable for 
# phylogenetic tree inference.
#
# The goal is to define a solid diploid backbone, which should be covered by
# outgroup sequences as well, and then use it to filter out polyploid 
# sequences with diploid block overlap < $MINBLOCKOVERLAP
#
# Edit variables below to define diploid, outgroup diploid and polyploid species.
#
# Only the longest sequence is taken for diploid species. 
#
# Note that outgroups are not used to compute diploid block.
#
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018

die "# usage: $0 <input MSA file> <output block MSA file>" if(!$ARGV[1]);

########################### user-defined variables ######################

my $VERBOSE = 0;

# min length of diploid aligned block
my $MINBLOCKLENGTH = 100;
my $MAXGAPSPERBLOCK = 100; # tolerated gaps for diploids in block
my $MINBLOCKOVERLAP = 0.50; # fraction of diploid block covered by outgroups and polyploid seqs
my $LONGTRINITYISOFPOLYPLOIDS = 1; # set to 1 if you wish that short trinity isoforms of the same polyploid sp 
                                   # are removed. Isoforms are identified with regex $ISOFORMREGEX
my $ISOFORMREGEX = 'c\d+_g\d+_i';

# sequences to be removed
my %to_be_removed = qw ( [syl_Cor_80_75] [syl_Gre_80_75] );

# user-defined hash of taxon names to be shortened
my %long2short = qw(
[Bdistachyon_314_v3.1] _Bdis
[Osativa_323_v7.0.transcript] _Osat
[Hvulgare_IBSC2016.HQ.cds] _Hvul
[arb_80_75] _Barb
[boi_80_75] _Bboi
[hyb_80_75] _Bhyb
[mex_80_75] _Bmex
[pho_80_75] _Bpho
[B422_80_75] _B422
[pin_80_75] _Bpin
[ret_80_75] _Bret
[rup_80_75] _Brup
[sta_80_75] _Bsta
[syl_Esp_80_75] _Bsyl
);

# short names of species used to define diploid block (see %long2short below)
# 1 per taxon will be selected for optimizing diploid block 
my @diploids = qw( _Bsta _Bdis _Barb _Bpin _Bsyl ); 

# diploid outgroups: best block-overlapping sequence will be conserved
my @outgroups = qw( _Osat _Hvul );

# polyploid species, only sequences ovarlapping the block will be conserved  
my @polyploids = qw( _Bhyb _Bboi _Bret _Bmex _Brup _Bpho _B422 );

# user-defined contribution of diploid species to block width calculations
# # arbitrary metric; in this case only one of Bpin or Bsyl is required
my %diploids4width = (
   '_Bsta'=>1,
   '_Bdis'=>1,
   '_Barb'=>1,
   '_Bpin'=>0.2,
   '_Bsyl'=>0.2
);

################################################################################

print "# params: MINBLOCKLENGTH=$MINBLOCKLENGTH MAXGAPSPERBLOCK=$MAXGAPSPERBLOCK\n";
print "# params: MINBLOCKOVERLAP=$MINBLOCKOVERLAP LONGTRINITYISOFPOLYPLOIDS=$LONGTRINITYISOFPOLYPLOIDS\n";
print "# diploids: ".join(',',@diploids)."\n";
print "# outgroups: ".join(',',@outgroups)."\n";
print "# polyploids: ".join(',',@polyploids)."\n\n";

################################################################################

my $MSAwidth = 0;
my $total_seqs_in_block = 0;
my $blockMSA = '';
my (%FASTA,%length,%inblock,%gaps,%header2taxon,%isoform);
my ($ngaps,$bases,$header,$seq,$previous_seq,$taxon,$short_taxon);
my ($coord,$coord3,$coord5,$overlap,$is_redundant);

my $input_MSA_file = $ARGV[0];
my $block_MSA_file = $ARGV[1];

## Simplify headers and exclude sequences to be removed 
open(MSA,"<",$input_MSA_file) || die "# ERROR: cannot open input file $input_MSA_file, exit\n";
while(my $line = <MSA>)
{
   if($line =~ /^>(.*?)(\[\S+?\])$/)
   {
		$header = $1;
      $taxon = $2;
		$short_taxon = '';
		
		if(!$to_be_removed{ $taxon })
		{ 
			$MSAwidth = 0;

			#shorten taxon string using user-defined %long2short
			$short_taxon = $long2short{$taxon} || "_$taxon";
		 	$header .= $short_taxon;
			$header2taxon{ $header } = $short_taxon;
		}
   }
   else
	{	
		if($short_taxon ne '')
		{	
			chomp($line);
			$FASTA{ $short_taxon }{ $header } .= $line;
			
			# record sequence length
         $bases = ($line =~ tr/[ACGTN]//);
			$length{ $header } += $bases;
			$MSAwidth += length($line);
		}
	}
}
close(MSA);

## define diploid block by checking the overlap of the longest sequences of taxons 
my $diploid_block_ok = 0;
my $diploid5prime = 0;
my $diploid3prime = $MSAwidth-1;
my $diploid_midcoord = 0;
my $diploid_block_length = 0;
foreach $taxon (@diploids)
{
	foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
	{	
		print "$taxon $length{$header} $header\n" if($VERBOSE);
		
		# to know which diploids have been considered for this calculation
		if($diploids4width{$taxon})
		{
			$diploid_block_ok += $diploids4width{$taxon};
		}

		# record 5' side of block of aligned diploids
      $coord = $diploid5prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord++ }
      if($coord > $diploid5prime){ $diploid5prime = $coord }

      # record 3' side of block of aligned diploids
      $coord = $diploid3prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord-- }
      if($coord < $diploid3prime){ $diploid3prime = $coord }
      print "# after $taxon : $diploid5prime to $diploid3prime\n" if($VERBOSE);

      $diploid_block_length = $diploid3prime - $diploid5prime + 1;

		$inblock{$header} = $total_seqs_in_block++;
		
		last; #consider only longest sequence of $short_taxon
	} 				
}

# compute block length
$diploid_block_length = $diploid3prime - $diploid5prime + 1;

if($diploid_block_ok <= 3)
{
   print "# MSA skipped as not all required diploids are aligned: $diploid_block_ok\n";
   exit;
}
elsif($diploid_block_length < $MINBLOCKLENGTH)
{
	print "# MSA skipped due to short diploid block: $diploid_block_length\n";
   exit;
}

# compute midpoint of diploid block
$diploid_midcoord = $diploid3prime - (($diploid3prime - $diploid5prime)/2);

print "# aligned diploid block: $diploid5prime < $diploid_midcoord > $diploid3prime\n\n";

## check gap fraction of sequences in diploid block
foreach $header (sort {$inblock{$a}<=>$inblock{$b}} keys(%inblock))
{
	$seq = substr( $FASTA{ $header2taxon{$header} } { $header },
							$diploid5prime,
                     $diploid3prime - $diploid5prime + 1 );

	$ngaps = ($seq =~ tr/\-//);

	if($ngaps > $MAXGAPSPERBLOCK)
	{
		print "# MSA skipped due to gappy diploid block: $ngaps ($header)\n";
		exit;
	}
}

## check overlap of outgroup sequences
foreach $taxon (@outgroups)
{
   foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
   {
      print "outg $taxon $length{$header} $header\n" if($VERBOSE);

		# make sure this sequence overlaps with diploid block
      # 5'
      $coord = $diploid5prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord++ }
      $coord5 = $coord;
      # 3'
      $coord = $diploid3prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord-- }
      $coord3 = $coord;

		$overlap = ($coord3 - $coord5 + 1) / $diploid_block_length;
	
		if( $overlap >= $MINBLOCKOVERLAP ){
			$inblock{ $header } = $total_seqs_in_block++;			
		}
		else
		{
			print "# MSA skipped as outgroup $taxon has poor block coverage: $overlap\n";
			exit;
		}

		last;
	}
}

## identify polyploid isoforms if requested
if($LONGTRINITYISOFPOLYPLOIDS)
{
	my ($isof);
	foreach $taxon (@polyploids)
	{
		my %tx_isoforms;
		foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
		{
			if($header =~ m/($ISOFORMREGEX)/)
			{
				$isoform = $1;
				$tx_isoforms{ $isoform }++;
				if($tx_isoforms{ $isoform } > 1)
				{
					$isoform{$taxon}{$header} = 1;
				}
			}
		}
	}
}

## check overlap of polyploid sequences
foreach $taxon (@polyploids)
{
	my (@seqs_in_block); # only for this poly species

   foreach $header (sort {$length{$b}<=>$length{$a}} keys(%{ $FASTA{$taxon} }))
   {
      print "poly $taxon $length{$header} $header\n" if($VERBOSE);
		
		if($isoform{$taxon}{$header})
		{
			print "# polyploid sequence $header skipped for being a redundant Trinity isoform\n";
			next;
		}

		# make sure this sequence overlaps with diploid block
      # 5'
      $coord = $diploid5prime;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord++ }
      $coord5 = $coord;
      # 3'
      $coord = $diploid3prime;;
      while(substr($FASTA{$taxon}{$header},$coord,1) eq '-'){ $coord-- }
      $coord3 = $coord;

		$overlap = ($coord3 - $coord5 + 1) / $diploid_block_length;

		if( $overlap >= $MINBLOCKOVERLAP )
		{
			print "$header $overlap\n" if($VERBOSE);

			# 1st seq of this poly species to be considered
			if(!@seqs_in_block)
			{
				$inblock{ $header } = $total_seqs_in_block++;
				push(@seqs_in_block, $header); # header of this sequence
			}
			else # make sure 2nd, 3rd seqs are not redundant
			{
				$is_redundant = 0;	
				$seq = substr( $FASTA{ $header2taxon{$header} } { $header },
                        $coord5,
                        $coord3 - $coord5 + 1 );				
				
				foreach my $previous_header (@seqs_in_block)
				{
					$previous_seq = substr( $FASTA{ $header2taxon{$header} } { $previous_header },
													$coord5,
													$coord3 - $coord5 + 1 );

					print "$seq $header\n$previous_seq $previous_header\n\n" if($VERBOSE);

					if($previous_seq eq $seq)
					{
						print "# polyploid redundant sequence skipped: $header\n";			
						$is_redundant = 1;
						last;
					}
				}	

				if($is_redundant == 0)
				{
					$inblock{ $header } = $total_seqs_in_block++;
					push(@seqs_in_block, $header);
				}
			}
      }
      else
      {
         print "# polyploid sequence $header skipped due to poor overlap: $overlap\n";
      }
	}
}


# create MSA of sequences trimmed to cover diploid block
open(BLOCKMSA,">",$block_MSA_file) || die "# ERROR: cannot create $block_MSA_file\n";

foreach $header (sort {$inblock{$a}<=>$inblock{$b}} keys(%inblock))
{
	print BLOCKMSA ">$header\n";
	print BLOCKMSA substr( $FASTA{ $header2taxon{$header} } { $header }, 
							$diploid5prime,
							$diploid3prime - $diploid5prime + 1 ) ."\n";		
}

close(BLOCKMSA);

print "\n# outfile: $block_MSA_file\n";
