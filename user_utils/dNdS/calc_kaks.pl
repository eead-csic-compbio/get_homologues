#!/usr/bin/perl
use strict;
use Statistics::Basic qw(:all nofill);

my $KAKSEXE = 'subopt-kaks-master/src/yn00_cds_prealigned';

my ($dir,$outdir) = (@ARGV);
if(!$ARGV[1]){ die "# need align and output dirs\n"; } 

opendir(DIR,$dir);
my @fnafiles = grep{/\.fna/} readdir(DIR);
closedir(DIR);

my ($n_of_aligns) = (0);

#header
print "<omega>\tpairs\tn_of_seqs\tfile\n";

foreach my $file (@fnafiles)
{			
   my $newname = $file;
	my $n_of_seqs = 0;
	open(FNA,$dir.$file);
	while(<FNA>)
	{
		if(/^>/)
		{
			$n_of_seqs++; 
		}
	}
	close(FNA);
	
   next if($n_of_seqs < 4);

   $newname =~ s/\.fna$/.tab/;

   #print "$KAKSEXE --output $outdir/$newname $dir/$file\n";
   system("$KAKSEXE --output $outdir/$newname $dir/$file") if(!-e "$outdir/$newname");

   my ($mean,$pairs,@omega) = ('n/a','n/a');
   open(TAB,"$outdir/$newname");
	while(<TAB>)
	{
		next if(/^SEQ/);
		my @data = split(/\t/,$_);

		next if($data[4] eq '-nan');

		if($data[4] == 99.00) 
		{		
			#next;
			$data[4] = 0.00; # convert pairs with dS==0 to omega = 0
		}

		push(@omega,$data[4]);
		$pairs++;
	}
	close(TAB);

	next if(scalar(@omega) < 1);

   $mean = mean(@omega);

	print "$mean\t$pairs\t$n_of_seqs\t$file\n";

	$n_of_aligns++;
}	

