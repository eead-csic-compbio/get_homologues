Authors
===============================================================================
GET_HOMOLOGUES is designed, created and maintained at the Laboratory of Computational 
Biology at Estacion Experimental de Aula Dei (EEAD/CSIC http://www.eead.csic.es/compbio) 
in Zaragoza, Spain, and at the Center for Genomic Sciences of Universidad Nacional 
Autonoma de Mexico (CCG/UNAM http://www.ccg.unam.mx/~vinuesa).

The program was written mostly by Bruno Contreras-Moreira and Pablo Vinuesa,
but includes code and binaries from other authors:

 OrthoMCL v1.4 (www.orthomcl.org , PubMed:12952885)
 mcl v14-137 (http://micans.org/mcl , PubMed=11917018)
 COGtriangles v2.1 (sourceforge.net/projects/cogtriangles , PubMed=20439257)
 NCBI Blast-2.2 (blast.ncbi.nlm.nih.gov , PubMed=9254694,20003500)
 BioPerl v1.5.2 (www.bioperl.org , PubMed=12368254)
 HMMER 3.1b2 (http://hmmer.org)
 Pfam (http://pfam.sanger.ac.uk , PubMed=19920124)
 PHYLIP 3.695 (http://evolution.genetics.washington.edu/phylip) 
 Transdecoder r20140704 (http://transdecoder.sf.net , PubMed=23845962)
 
License
===============================================================================
Please check the LICENSE file.

Key literature references:
================================================================================
1) Contreras-Moreira,B. and Vinuesa,P. (2013) GET_HOMOLOGUES, a versatile software
package for scalable and robust microbial pangenome analysis. Appl.Environ.Microbiol. 
79:7696-7701.

2) Vinuesa P and Contreras-Moreira B (2015) Robust Identification of Orthologues 
and Paralogues for Microbial Pan-Genomics Using GET_HOMOLOGUES: A Case Study of 
pIncA/C Plasmids. In Bacterial Pangenomics, Methods in Molecular Biology Volume 
1231, 203-232, edited by A Mengoni, M Galardini and M Fondi.

3) Li L, Stoeckert CJ Jr, Roos DS (2003) OrthoMCL: identification of ortholog 
groups for eukaryotic genomes. Genome Res. 13(9):2178-89.

4) Kristensen DM, Kannan L, Coleman MK, Wolf YI, Sorokin A, Koonin EV, 
Mushegian A (2010) A low-polynomial algorithm for assembling clusters of orthologous 
groups from intergenomic symmetric best matches. Bioinformatics 26(12):1481-7.

5) Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W and Lipman DJ 
(1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search 
programs. Nucl. Acids Res. 25(17): 3389-3402. 

6) Stajich JE, Block D, Boulez K, Brenner SE, Chervitz SA, Dagdigian C, Fuellen G, 
Gilbert JG, Korf I, Lapp H, Lehväslaiho H, Matsalla C, Mungall CJ, Osborne BI, 
Pocock MR, Schattner P, Senger M, Stein LD, Stupka E, Wilkinson MD, Birney E. (2002)
The Bioperl toolkit: Perl modules for the life sciences. Genome Res. 12(10):1611-8.

7) hmmscan :: search sequence(s) against a profile database HMMER 3.0 (March 2010) 
http://hmmer.org Copyright (C) 2010 Howard Hughes Medical Institute.
Freely distributed under the GNU General Public License (GPLv3).

8) Finn RD, Bateman A, Clements J, Coggill P, Eberhardt RY, Eddy SR, Heger A, 
Hetherington K, Holm L, Mistry J, Sonnhammer EL, Tate J, Punta M. (2014)
Pfam: the protein families database. Nucleic Acids Res. 42:D222-30. 

9) Felsenstein, J. 2005. PHYLIP (Phylogeny Inference Package) version 3.6. Distributed 
by the author. Department of Genome Sciences, University of Washington, Seattle.

10) Haas BJ, Papanicolaou A, Yassour M et al. (2013) De novo transcript sequence 
reconstruction from RNA-seq using the Trinity platform for reference generation 
and analysis. Nat Protoc. 8(8):1494-512.

System requirements and installation:
================================================================================
This software has been tested in Linux and MacOSX systems. It includes a 
a Perl script and compiled binaries from the abovementioned authors.
In order to use this software follow these steps:

1) Unpack the GET_HOMOLOGUES software package with 'tar xvfz get_homologues_Y.X.tgz' 

2) cd get_homologues_Y.X/

3) perl install.pl

4) ./get_homologues.pl -v 

5) ./get_homologues.pl -h OR ./get_homologues-est.pl -h

Further documentation 
================================================================================
Please check the accompanying manuals in PDF format for more details.

