#!/usr/bin/perl
use strict;    

my @barleys = qw( 73 SC Hs );

my %annotfiles = ( 
'SBCC073_accessory.pfam.enrich.tab'=>'73',
'Scarlett_accessory.pfam.enrich.tab'=>'SC',
'spontaneum_accessory.pfam.enrich.tab'=>'Hs'
);

my (%stats,%adjP,%PfamDescr);
my ($Pfam,$exp,$ctrl,$dum,$dum,$dum,$adjP,$descr,$abbrev);

foreach my $file (keys(%annotfiles))
{ 
  $abbrev = $annotfiles{$file};
  #print "# parsing $abbrev $file\n";

  open(ANNOT,$file);
  while(<ANNOT>)
  {
    ##PfamID  counts(exp) counts(ctr) freq(exp) freq(ctr) p-value p-value(adj)  description
    #PF00085 7 61  1.913e-02 2.324e-03 4.594e-05 5.972e-02 Thioredoxin
    next if(/^#/);
    chomp;
    ($Pfam,$exp,$ctrl,$dum,$dum,$dum,$adjP,$descr) = split(/\t/,$_);

    if(!$PfamDescr{$Pfam}){ $PfamDescr{$Pfam} = $descr }

    $stats{$Pfam}{'total'} += $exp;
    $stats{$Pfam}{$abbrev} = $exp;
    $adjP{$Pfam}{$abbrev} = $adjP;
  }
}
close(ANNOT);


print "#".join("\t",@barleys)."\tadjP\tadjP\tadjP\tdomain\n";

foreach $Pfam (sort {$stats{$b}{'total'}<=>$stats{$a}{'total'}} keys(%stats))
{	
  #print "$Pfam\t$stats{$Pfam}{'total'}";
  foreach $abbrev (@barleys)
  {
    printf("%d\t",$stats{$Pfam}{$abbrev} || 0); 
  }
  foreach $abbrev (@barleys)
  {
    printf("%s\t",$adjP{$Pfam}{$abbrev} || 'NA');
  }
  print "$Pfam $PfamDescr{$Pfam}\n";
}

#data = read.table(file="accessory_stats_edited.tab",header=T,sep="\t")
#rnames <- data[,7]
#mat_data <- data.matrix(data[,1:3])
#rownames(mat_data) <- rnames
#colnames(mat_data) <- c("73","SC","Hs")
#palette <- colorRampPalette(c("lightyellow","red"))(n = 50)
#pdf("accessory_stats_edited.pdf")
#heatmap.2(mat_data,cellnote=mat_data,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",dendrogram="none",margins =c(1,25),colsep=c(0:3),rowsep=c(0:50),sepwidth=c(0.01,0.01), sepcolor="grey", col=palette, lhei = c(0.05,5),labCol='')
#dev.off()

