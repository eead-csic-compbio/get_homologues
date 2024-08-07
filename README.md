## GET_HOMOLOGUES 

A versatile software package for pan-genome analysis, including GET_HOMOLOGUES and GET_HOMOLOGUES-EST. It includes several algorithms designed for:
 * Clustering coding sequences in homologous (possibly orthologous) groups, on the grounds of sequence similarity. By default GET_HOMOLOGUES compares protein sequences, while GET_HOMOLOGUES-EST aligns nucleotide sequences which can be CDS or transcripts.
 * Definition of pan- and core-genomes by calculation of overlapping sets of proteins or CDSs.

GET_HOMOLOGUES has been used mostly with bacterial data (see citing 
[papers](https://scholar.google.es/scholar?start=0&hl=en&as_sdt=2005&cites=5259912818944685430)).
 
Instead, GET_HOMOLOGUES-EST has been used mostly with plants (see citing
[papers](https://scholar.google.es/scholar?oi=bibs&hl=en&cites=14330917787074873427&as_sdt=5) and 
was originally benchmarked with genomes and transcriptomes of *Arabidopsis thaliana* and *Hordeum vulgare*, 
available at [http://floresta.eead.csic.es/plant-pan-genomes](http://floresta.eead.csic.es/plant-pan-genomes), 
and used to analyze the pan-genomes of [*Brachypodium distachyon*](https://brachypan.jgi.doe.gov) and 
[*Brachypodium hybridum*](http://floresta.eead.csic.es/plant-pan-genomes/Bhybridum) ([press release](https://jgi.doe.gov/more-the-merrier-making-case-for-plant-pan-genomes)).

- [Installation](#installation)
- [Documentation](#documentation)
- [Citation](#citation)
- [Credits](#rredits)
- [Graphical summary](#graphical-summary)
- [Related software](#related-software)
- [Bugs](#bugs)
- [Funding](#funding)

## Installation

[![Build Status](https://app.travis-ci.com/eead-csic-compbio/get_homologues.svg?branch=master)](https://app.travis-ci.com/eead-csic-compbio/get_homologues)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/get_homologues/badges/version.svg)](https://anaconda.org/bioconda/get_homologues)
[![DockerHub](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/csicunam/get_homologues)

Installation instructions, including the [bioconda package](https://anaconda.org/bioconda/get_homologues), are available in the
[manual](http://eead-csic-compbio.github.io/get_homologues/manual/manual.html#SECTION00030000000000000000)
and the [README.txt](./README.txt) file.

Check also the [Docker image](https://hub.docker.com/r/csicunam/get_homologues).

## Documentation

Manuals are available at:

|version|HTML|
|-------|----|
|original, for the analysis of bacterial pan-genomes|[manual](http://eead-csic-compbio.github.io/get_homologues/manual/)|
|EST, for the analysis of intra-species eukaryotic pan-genomes, tested on plants|[manual-est](http://eead-csic-compbio.github.io/get_homologues/manual-est/)|

In addition, there are two **tutorials** are available:

* [Pangenome analysis of plant transcripts and coding sequences](http://eead-csic-compbio.github.io/get_homologues/plant_pangenome/protocol.html), published in 2022.

* [From genomes to pangenomes: understanding variation among individuals and species](http://eead-csic-compbio.github.io/get_homologues/tutorial/pangenome_tutorial.html), which includes step by step instructions for both bacterial and plant data, first released in 2017.

## Citation

The original GET_HOMOLOGUES, **suitable for bacterial genomes**, was described in:

[Contreras-Moreira B, Vinuesa P (2013) Appl. Environ. Microbiol. 79:7696-7701](http://aem.asm.org/content/79/24/7696.long)

[Vinuesa P, Contreras-Moreira B (2015) Methods in Molecular Biology Volume 1231, 203-232](http://link.springer.com/protocol/10.1007%2F978-1-4939-1720-4_14)

GET_HOMOLOGUES-EST, adapted to the study of **intra-specific eukaryotic pan-genomes and pan-transcriptomes**, was described in:

[Contreras-Moreira B, Cantalapiedra CP et al (2017) Front. Plant Sci. 10.3389/fpls.2017.00184](http://journal.frontiersin.org/article/10.3389/fpls.2017.00184/full)

[Contreras-Moreira B, Rodriguez del Rio A et al (2022) Methods in Molecular Biology https://doi.org/10.1007/978-1-0716-2429-6_9](https://link.springer.com/protocol/10.1007/978-1-0716-2429-6_9)

## Credits 

GET_HOMOLOGUES is designed, created and maintained at the Computational and 
Structural Biology group at Estacion Experimental de Aula Dei (EEAD-CSIC 
https://www.eead.csic.es/compbio) in Zaragoza, Spain, and at the Center for 
Genomic Sciences of Universidad Nacional Autonoma de Mexico (CCG/UNAM 
http://www.ccg.unam.mx/~vinuesa).

The program was written mostly by Bruno Contreras-Moreira and Pablo Vinuesa,
with contributions from Carlos P Cantalapiedra, Alvaro Rodriguez del Rio, Ruben
Sancho, Roland Wilhelm, David A Wilkinson and many others (see [CHANGES.txt](./CHANGES.txt)). 
It also includes code and binaries from other authors:

 OrthoMCL v1.4 (www.orthomcl.org , PubMed:12952885)
 mcl v14-137 (http://micans.org/mcl , PubMed=11917018)
 COGtriangles v2.1 (sourceforge.net/projects/cogtriangles , PubMed=20439257)
 NCBI Blast-2.16.1+ (blast.ncbi.nlm.nih.gov , PubMed=9254694,20003500)
 BioPerl v1.5.2 (www.bioperl.org , PubMed=12368254)
 HMMER 3.1b2 (http://hmmer.org)
 Pfam (http://pfam.xfam.org , PubMed=19920124)
 PHYLIP 3.695 (http://evolution.genetics.washington.edu/phylip) 
 Transdecoder r20140704 (http://transdecoder.github.io , PubMed=23845962)
 MVIEW 1.60.1 (https://github.com/desmid/mview, PubMed=9632837)
 diamond 0.8.25 (https://github.com/bbuchfink/diamond, PubMed=25402007)

## Graphical summary

![**Legend.** Main features of GET_HOMOLOGUES.](./pics/summary.jpg)

![**Legend.** Flowchart and features of GET_HOMOLOGUES-EST.](./pics/EST.jpg)

## Related software

[GET_PHYLOMARKERS](https://github.com/vinuesa/get_phylomarkers) uses 
twin nucleotide & peptide clusters produced by GET_HOMOLOGUES to compute robust multi-gene and pangenome phylogenies.
Check the [manual](https://vinuesa.github.io/get_phylomarkers), the [tutorial](https://link.springer.com/protocol/10.1007/978-1-0716-2429-6_9),
and the [Docker image](https://hub.docker.com/r/csicunam/get_homologues).

A related piece of software was released in 2023 called [GET_PANGENES](https://github.com/Ensembl/plant-scripts/tree/master/pangenes),
which takes FASTA and GFF files as input and explicitely considers gene collinearity.

## Bugs

The code is regularly patched (see [CHANGES.txt](./CHANGES.txt) in each release. We kindly ask you to report errors or bugs as 
[GitHub issues](https://github.com/eead-csic-compbio/get_homologues/issues)
and to acknowledge the use of the software in scientific publications.

## Funding

[GET_HOMOLOGUES](https://bio.tools/get_homologues) is part of the [INB/ELIXIR-ES](https://inb-elixir.es) resources portfolio:

![logo_ELIXIRES](pics/logoELIXIRES.png)

*Funding:* Fundacion ARAID, Consejo Superior de Investigaciones Cientificas, DGAPA-PAPIIT UNAM, CONACyT, FEDER, MINECO, DGA-Obra Social La Caixa.

![logo CSIC](pics/logoCSIC.png) ![logo ARAID](pics/logoARAID.gif) ![logo UNAM](pics/logoUNAM.png)

