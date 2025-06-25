Authors
===============================================================================
GET_HOMOLOGUES is designed, created and maintained at the Computational and 
Structural Biology group at Estacion Experimental de Aula Dei (EEAD-CSIC 
https://www.eead.csic.es/compbio) in Zaragoza, Spain, and at the Center for 
Genomic Sciences of Universidad Nacional Autonoma de Mexico (CCG/UNAM 
http://www.ccg.unam.mx/~vinuesa).

The program was written mostly by Bruno Contreras-Moreira and Pablo Vinuesa,
with contributions from Carlos P Cantalapiedra, Alvaro Rodriguez del Rio, Ruben
Sancho, Roland Wilhelm, David A Wilkinson and many others (see CHANGES.txt). 
It also includes code and binaries from other authors:

 OrthoMCL v1.4 (www.orthomcl.org , PubMed:12952885)
 mcl v14-137 (http://micans.org/mcl , PubMed=11917018)
 COGtriangles v2.1 (sourceforge.net/projects/cogtriangles , PubMed=20439257)
 NCBI Blast-2.16.0+ (blast.ncbi.nlm.nih.gov , PubMed=9254694,20003500)
 BioPerl v1.5.2 (www.bioperl.org , PubMed=12368254)
 HMMER 3.1b2 (http://hmmer.org)
 Pfam (http://pfam.xfam.org , PubMed=19920124)
 PHYLIP 3.695 (https://phylipweb.github.io/phylip) 
 Transdecoder r20140704 (http://transdecoder.github.io , PubMed=23845962)
 MVIEW 1.60.1 (https://github.com/desmid/mview, PubMed=9632837)
 diamond 0.8.25 (https://github.com/bbuchfink/diamond, PubMed=25402007)
 
License
===============================================================================
Please check the LICENSE file.

Key literature references:
================================================================================

1) Contreras-Moreira B, Cantalapiedra CP, Garcia Pereira MJ, Gordon S, Vogel JP, 
Igartua E, Casas AM and Vinuesa P (2017) Analysis of plant pan-genomes and 
transcriptomes with GET_HOMOLOGUES-EST, a clustering solution for sequences of 
the same species. Front. Plant Sci. 10.3389/fpls.2017.00184

2) Contreras-Moreira B and Vinuesa P (2013) GET_HOMOLOGUES, a versatile software
package for scalable and robust microbial pangenome analysis. Appl.Environ.Microbiol. 
79:7696-7701.

3) Contreras-Moreira B, Rodriguez del Rio A, Cantalapiedra CP, Sancho R, Vinuesa P
(2022) Pangenome Analysis of Plant Transcripts and Coding Sequences. In: Pereira-Santana A, 
Gamboa-Tuz SD, Rodriguez-Zapata LC (eds) Plant Comparative Genomics. Methods 
in Molecular Biology, vol 2512. Humana, New York, NY.  

4) Vinuesa P and Contreras-Moreira B (2015) Robust Identification of Orthologues 
and Paralogues for Microbial Pan-Genomics Using GET_HOMOLOGUES: A Case Study of 
pIncA/C Plasmids. In Bacterial Pangenomics, Methods in Molecular Biology Volume 
1231, 203-232, edited by A Mengoni, M Galardini and M Fondi.

5) Li L, Stoeckert CJ Jr, Roos DS (2003) OrthoMCL: identification of ortholog 
groups for eukaryotic genomes. Genome Res. 13(9):2178-89.

6) Kristensen DM, Kannan L, Coleman MK, Wolf YI, Sorokin A, Koonin EV, 
Mushegian A (2010) A low-polynomial algorithm for assembling clusters of orthologous 
groups from intergenomic symmetric best matches. Bioinformatics 26(12):1481-7.

7) Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W and Lipman DJ 
(1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search 
programs. Nucl. Acids Res. 25(17): 3389-3402. 

8) Stajich JE, Block D, Boulez K, Brenner SE, Chervitz SA, Dagdigian C, Fuellen G, 
Gilbert JG, Korf I, Lapp H, Lehvaslaiho H, Matsalla C, Mungall CJ, Osborne BI, 
Pocock MR, Schattner P, Senger M, Stein LD, Stupka E, Wilkinson MD, Birney E. (2002)
The Bioperl toolkit: Perl modules for the life sciences. Genome Res. 12(10):1611-8.

9) hmmscan :: search sequence(s) against a profile database HMMER 3.0 (March 2010) 
http://hmmer.org Copyright (C) 2010 Howard Hughes Medical Institute.
Freely distributed under the GNU General Public License (GPLv3).

10) Finn RD, Bateman A, Clements J, Coggill P, Eberhardt RY, Eddy SR, Heger A, 
Hetherington K, Holm L, Mistry J, Sonnhammer EL, Tate J, Punta M. (2014)
Pfam: the protein families database. Nucleic Acids Res. 42:D222-30. 

11) Felsenstein, J. 2005. PHYLIP (Phylogeny Inference Package) version 3.6. Distributed 
by the author. Department of Genome Sciences, University of Washington, Seattle.

12) Haas BJ, Papanicolaou A, Yassour M et al. (2013) De novo transcript sequence 
reconstruction from RNA-seq using the Trinity platform for reference generation 
and analysis. Nat Protoc. 8(8):1494-512.

13) Brown NP, Leroy C, Sander C (1998) MView: A Web compatible database search or 
multiple alignment viewer. Bioinformatics. 14 (4):380-381. 

14) Buchfink B, Xie C, Huson DH (2015) Fast and sensitive protein alignment using 
DIAMOND. Nat Methods. 12(1):59-60

System requirements and installation:
================================================================================
This software has been tested in Linux and MacOSX systems. It includes a 
a Perl script and compiled binaries from the abovementioned authors.

In order to install and test this software please follow these steps:

 Download a bundled release from https://github.com/eead-csic-compbio/get_homologues/releases

  $ tar xvfz get_homologues_X.Y.tgz # unpack it
  $ cd get_homologues_X.Y
  $ perl install.pl # follow the indications in case some required part is missing.

  $ ./get_homologues.pl -v # will tell exactly which features are available.
  $ ./get_homologues.pl -d sample_buch_fasta # test, output should be like sample_output.txt.

  Optionally modify your $PATH environment variable to include get_homologues.pl. 
  Please copy the following lines to the .bash_profile or .bashrc files, found in your home directory, 
  replacing [INSTALL_PATH] by the full path of the installation folder:

  $ export GETHOMS=[INSTALL_PATH]/get_homologues_X.Y
  $ export PATH=${GETHOMS}/:${PATH}

  This change will be effective in a new terminal or after running: 
  $ source ~/.bash_profile

 If you prefer a copy of the software that can be updated in the future you can get it from the GitHub repository with:

  $ git clone https://github.com/eead-csic-compbio/get_homologues.git
  $ perl install.pl
 
 You would then be able to update it at anytime with:
  $ cd get_homologues
  $ git pull

 Finally, you can also install the software from bioconda as follows:

  $ conda activate bioconda
  $ conda create -n get_homologues -c conda-forge -c bioconda get_homologues
  $ conda activate get_homologues

  # only if you want to install Pfam or SwissProt db
  $ perl install.pl


Further documentation 
================================================================================
Please check the accompanying manuals in PDF format for more details, or the 
online tutorials and manuals at https://github.com/eead-csic-compbio/get_homologues

