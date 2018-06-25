norm_bitscore.pl
================

 This Perl5 script takes a file with tab-separated BLASTP/BLASTN results, 
with the custom-format used by get_homologues, and computes length-normalized 
E-value & bitscores, as done originally in software OrthoFinder 
(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2)

Authors: Alvaro Rodriguez del Rio and Bruno Contreras-Moreira 

## Example

You can test this script with the provided sample file as follows:

```
perl norm_bitscore.pl -i sample.blastp.gz -o sample.blastp.norm
```

## Comparing raw and normalized BLAST scores with your data

1) Calculate clusters with raw scores:

```
perl get_homologues-est.pl -d your_sequences -M -t 0 &> log.OMCL.raw
```

2) Let's say this created a folder with cluster named *your_sequences_homologues/clusters_0taxa_algOMCL_e0_*.
We'll renamed them so that they're not overwritten:

```
mv your_sequences_homologues/clusters_0taxa_algOMCL_e0_ \
  your_sequences_homologues/clusters_0taxa_algOMCL_e0_raw
mv your_sequences_homologues/clusters_0taxa_algOMCL_e0_.cluster_list \
  your_sequences_homologues/clusters_0taxa_algOMCL_e0_raw.cluster_list
```

3) Save raw BLAST outfiles and put them away in a separate folder named *raw/*:

```
cd your_sequences_homologues/
mkdir raw norm
cp *.blast.gz raw/
```

4) Compute normalized BLAST scores and put them away in a separate folder named *norm/*:

```
for FILE in *.blast.gz; do perl norm_bitscore.pl -i ori/$FILE -o norm/${FILE%.*}; done
gzip norm/*
rm -f *.blast.gz tmp/all.b*
ln -s norm/* .
cd tmp
mkdir raw
mv orthologues_010_* inparalogues_* raw/
cd ../..
```

5) Calculate clusters with normalized scores and compare them to the previous ones:

```
perl get_homologues-est.pl -d your_sequences -m cluster -M -t 0 &> log.OMCL.norm

mv your_sequences_homologues/clusters_0taxa_algOMCL_e0_ \
  your_sequences_homologues/clusters_0taxa_algOMCL_e0_norm
mv your_sequences_homologues/clusters_0taxa_algOMCL_e0_.cluster_list  
  your_sequences_homologues/clusters_0taxa_algOMCL_e0_norm.cluster_list

perl compare_clusters.pl -m -o intersection -d \
  your_sequences_homologues/clusters_0taxa_algOMCL_e0_raw,your_sequences_homologues/clusters_0taxa_algOMCL_e0_norm \
  &> log.intersection
```
