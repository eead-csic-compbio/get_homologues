## GET_HOMOLOGUES: a versatile software package for pan-genome analysis

[![Build Status](https://app.travis-ci.com/eead-csic-compbio/get_homologues.svg?branch=master)](https://app.travis-ci.com/eead-csic-compbio/get_homologues)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/get_homologues/badges/version.svg)](https://anaconda.org/bioconda/get_homologues)
[![DockerHub](https://badgen.net/badge/icon/docker?icon=docker&label)](https://hub.docker.com/r/csicunam/get_homologues)


This software is maintained by Bruno Contreras-Moreira (bcontreras _at_ eead.csic.es) and Pablo Vinuesa (vinuesa _at_ ccg.unam.mx). 
The original version, **suitable for bacterial genomes**, was described in:

![**Legend.** Main features of GET_HOMOLOGUES.](./pics/summary.jpg)

[Contreras-Moreira B, Vinuesa P (2013) Appl. Environ. Microbiol. 79:7696-7701](http://aem.asm.org/content/79/24/7696.long)

[Vinuesa P, Contreras-Moreira B (2015) Methods in Molecular Biology Volume 1231, 203-232](http://link.springer.com/protocol/10.1007%2F978-1-4939-1720-4_14)

The software was subsequently adapted to the study of **intra-specific eukaryotic pan-genomes** resulting in script GET_HOMOLOGUES-EST, described in:

![**Legend.** Flowchart and features of GET_HOMOLOGUES-EST.](./pics/EST.jpg)

[Contreras-Moreira B, Cantalapiedra CP et al (2017) Front. Plant Sci. 10.3389/fpls.2017.00184](http://journal.frontiersin.org/article/10.3389/fpls.2017.00184/full)

[Contreras-Moreira B, Rodriguez del Rio A et al (2022) Methods in Molecular Biology https://doi.org/10.1007/978-1-0716-2429-6_9](https://link.springer.com/protocol/10.1007/978-1-0716-2429-6_9)

GET_HOMOLOGUES-EST was benchmarked with genomes and transcriptomes of *Arabidopsis thaliana* and *Hordeum vulgare*, available at [http://floresta.eead.csic.es/plant-pan-genomes](http://floresta.eead.csic.es/plant-pan-genomes), and used to analyze the pan-genomes of [*Brachypodium distachyon*](https://brachypan.jgi.doe.gov) and [*Brachypodium hybridum*](http://floresta.eead.csic.es/plant-pan-genomes/Bhybridum) ([press release](https://jgi.doe.gov/more-the-merrier-making-case-for-plant-pan-genomes)).


Two **tutorials** are available:
* [Pangenome analysis of plant transcripts and coding sequences](http://eead-csic-compbio.github.io/get_homologues/plant_pangenome/protocol.html), published in 2022.

* [From genomes to pangenomes: understanding variation among individuals and species](http://eead-csic-compbio.github.io/get_homologues/tutorial/pangenome_tutorial.html), which includes step by step instructions for both bacterial and plant data, first released in 2017.
 

**Installation instructions**, including the [bioconda package](https://anaconda.org/bioconda/get_homologues), are available in the 
[manual](http://eead-csic-compbio.github.io/get_homologues/manual/manual.html#SECTION00030000000000000000) 
and the [README.txt](./README.txt) file.

|version|HTML|
|-------|----|
|original, for the analysis of bacterial pan-genomes|[manual](http://eead-csic-compbio.github.io/get_homologues/manual/)|
|EST, for the analysis of intra-species eukaryotic pan-genomes, tested on plants|[manual-est](http://eead-csic-compbio.github.io/get_homologues/manual-est/)|

A **[Docker image](https://hub.docker.com/r/csicunam/get_homologues)** is also available with GET_HOMOLOGUES 
bundled with [GET_PHYLOMARKERS](https://github.com/vinuesa/get_phylomarkers), ready to use. 
The GET_PHYLOMARKERS [manual](https://vinuesa.github.io/get_phylomarkers) 
explains how to use nucleotide & peptide clusters produced by GET_HOMOLOGUES to compute robust multi-gene and pangenome phylogenies.

A related piece of software was released in 2023 called [GET_PANGENES](https://github.com/Ensembl/plant-scripts/tree/master/pangenes),
which takes FASTA and GFF files as input and explicitely considers gene collinearity.

The code is regularly patched (see [CHANGES.txt](./CHANGES.txt) in each release, and has been used in a variety of studies 
(see citing papers [here](https://scholar.google.es/scholar?start=0&hl=en&as_sdt=2005&cites=5259912818944685430) and 
[here](https://scholar.google.es/scholar?oi=bibs&hl=en&cites=14330917787074873427&as_sdt=5), respectively).

We kindly ask you to report errors or bugs to the authors and to acknowledge the use of the software in scientific publications.

[GET_HOMOLOGUES](https://bio.tools/get_homologues) is part of the [INB/ELIXIR-ES](https://inb-elixir.es) resources portfolio:

![logo_ELIXIRES](pics/logoELIXIRES.png)

*Funding:* Fundacion ARAID, Consejo Superior de Investigaciones Cientificas, DGAPA-PAPIIT UNAM, CONACyT, FEDER, MINECO, DGA-Obra Social La Caixa.

![logo CSIC](pics/logoCSIC.png) ![logo ARAID](pics/logoARAID.gif) ![logo UNAM](pics/logoUNAM.png)

