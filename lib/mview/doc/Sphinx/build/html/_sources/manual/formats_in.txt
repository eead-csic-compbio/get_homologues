.. _ref_input_formats:

===============
 Input formats
===============

The code has been tested for the following formats and versions for protein
and nucleotide sequences:


NCBI BLAST family
=================

Sequence database search programs. See :ref:`ref_blast_rules` for details of
`BLAST` processing by MView.


BLAST+
------

Modern BLAST can generate output in a number of formats; run `blastp -help` to
see what is available.

MView can process two of these: the default BLAST output format (also obtained
with the BLAST command line option `-outfmt 0`) and a commented tabular format
produced with `-outfmt 7`. MView support and testing is as follows.

MView option: ``-in blast`` (with the default BLAST format)

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `blastp`     2.2.25+, 2.2.28+, 2.2.31+              ok
 `blastn`     2.2.25+, 2.2.28+, 2.2.31+              ok
 `blastx`     2.2.28+, 2.2.31+                       ok
 `tblastn`    2.2.28+, 2.2.31+                       ok
 `tblastx`    2.2.28+, 2.2.31+                       ok
 `psiblast`   2.2.28+, 2.2.31+                       ok
 ===========  ====================================== ============

MView option: ``-in blast`` (with `-outfmt 7` tabular BLAST format)

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `blastp`     2.2.28+, 2.2.31+                       ok
 `blastn`     2.2.28+, 2.2.31+                       ok
 `blastx`     2.2.28+, 2.2.31+                       ok
 `tblastn`    2.2.28+, 2.2.31+                       ok
 `tblastx`    2.2.28+, 2.2.31+                       ok
 `psiblast`   2.2.28+, 2.2.31+                       ok
 ===========  ====================================== ============

The `-outfmt 7` fields or columns recognised by MView are:

- the `std` defaults (`qseqid`, `sseqid`, `pident`, `length`, `mismatch`,
  `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`);
- `qseq`, `sseq` to report the query and sbjct sequences;
- `stitle` or `salltitles` to report the sbjct descriptions;
- `qframe` to report the query frame;
- `sframe` to report the sbjct frame.

Suitable strings of field specifiers for input to MView are:

    `blastp   -outfmt '7 std qseq sseq stitle'`

    `blastn   -outfmt '7 std qseq sseq stitle'`

    `blastx   -outfmt '7 std qseq sseq stitle qframe'`

    `tblastn  -outfmt '7 std qseq sseq stitle sframe'`

    `tblastx  -outfmt '7 std qseq sseq stitle qframe sframe'`

    `psiblast -outfmt '7 std qseq sseq stitle'`

MView will silently ignore any fields other than those listed above. The field
ordering given does not matter as MView identifies this automatically.

Note: There are some differences in the output (seen in BLAST `2.2.28+`)
produced by the default BLAST format versus `-outfmt 7`. The latter lacks the
list of ranked database hits and best scores and their associated fragment
counts, N. Unfortunately this can obscure interpretation of the results when
there are many alternative alignments for a given database hit (see
:ref:`ref_blast_rules`) so that in certain cases the MView output will be
slightly different; this is most noticeable with `tblastx`.

Note: A strange behaviour has been seen with `psiblast 2.2.28+`, where the
first database search produces many more hits than subsequent iterations
including the final iteration when run using the default output format. Using
`-outfmt 7` on the exact same search produces more consistent sets of hits.


BLAST series 2.2
----------------

MView option: ``-in blast``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `blastp`     2.2.5, 2.2.6                           ok
 `blastn`     2.2.6                                  ok
 `blastx`     2.2.6                                  ok
 `tblastn`    2.2.6                                  ok
 `tblastx`    2.2.6                                  ok
 `psiblast`   2.2.6                                  ok
 `phiblast`   2.2.6                                  ok
 ===========  ====================================== ============


BLAST series 2.0
----------------

MView option: ``-in blast``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `blastp`     2.0.4, 2.0.5, 2.0.9, 2.0.10            ok
 `blastn`     2.0.4, 2.0.5, 2.0.9, 2.0.14            ok
 `blastx`     2.0.5, 2.0.9                           ok
 `tblastn`    2.0.5, 2.0.10                          ok
 `tblastx`    2.0.5                                  ok
 `psiblast`   2.0.2, 2.0.4, 2.0.5, 2.0.6, 2.0.10     ok
 `phiblast`   2.0.9                                  ok
 ===========  ====================================== ============


BLAST series 1.4
----------------

MView option: ``-in blast``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `blastp`     1.4.7, 1.4.9                           ok
 `blastn`     1.4.9                                  ok
 `blastx`     1.4.9                                  ok
 `tblastn`    1.4.9                                  ok
 `tblastx`    1.4.9                                  ok
 ===========  ====================================== ============


------------------------------------------------------------------------------

WashU-BLAST family
==================

Sequence database search programs.


WU-BLAST series 2.0
-------------------

MView option: ``-in blast``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `blastp`     2.0a13, 2.0a19, 2.0                    ok
 `blastn`     2.0a19, 2.0                            ok
 `blastx`     2.0a19, 2.0                            ok
 `tblastn`    2.0a19, 2.0                            ok
 `tblastx`    2.0a19, 2.0                            ok
 ===========  ====================================== ============


------------------------------------------------------------------------------

Uni. Virginia FASTA family
==========================

Sequence database search programs:


FASTA series 36
---------------

MView option: ``-in uvfasta``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `fasta`      36.07, 36.3.3, 36.3.6, 36.3.7b         ok
 `fastx`      36.3.3, 36.3.6, 36.3.7b                ok
 `fasty`      36.3.6, 36.3.7b                        ok
 `tfastx`     36.3.3, 36.3.6, 36.3.7b                ok
 `tfasty`     36.3.3, 36.3.6, 36.3.7b                ok
 `ssearch`    36.3.3, 36.3.6, 36.3.7b                ok
 `ggsearch`   36.3.3, 36.3.6, 36.3.7b                ok
 `glsearch`   36.3.3, 36.3.6, 36.3.7b                ok
 `fastm`      35.3.5a, 36.3.6, 36.3.7b               experimental
 `fasts`      36.3.6, 36.3.7b                        experimental
 `fastf`      36.3.6, 36.3.7b                        experimental
 `tfastm`     36.3.7b                                experimental
 `tfasts`     36.3.7b                                experimental
 `tfastf`     36.3.7b                                experimental
 ===========  ====================================== ============


FASTA series 35
---------------

MView option: ``-in uvfasta``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `fasta`      35.04                                  ok
 `tfastx`     35.04                                  ok
 `ssearch`    35.04                                  ok
 `ggsearch`   35.04                                  ok
 `glsearch`   35.04                                  ok
 ===========  ====================================== ============


FASTA series 34
---------------

MView option: ``-in uvfasta``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `fasta34`    34.26.5                                ok
 ===========  ====================================== ============


FASTA series 3.0 - 3.4
----------------------

MView option: ``-in uvfasta``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `fasta`      3.0t76, 3.1t07, 3.2t01, 3.2t05, 3.2t07
 \            3.2t05, 3.2t07, 3.3t01, 3.3t07, 3.4t23 ok
 `fastx`      3.4t23                                 ok
 `fasty`      3.4t23                                 ok
 `tfasta`     3.4t23                                 ok
 `tfastx`     3.0t82, 3.4t23                         ok
 `tfasty`     3.4t23                                 ok
 `tfastxy`    3.1t07                                 ok
 ===========  ====================================== ============


FASTA series 2
--------------

MView option: ``-in uvfasta``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `fasta`      2.0u                                   ok
 `tfastx`     2.0u63                                 ok
 ===========  ====================================== ============


FASTA series 1
--------------

MView option: ``-in uvfasta``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `fasta`      1.6c24                                 ok
 ===========  ====================================== ============


------------------------------------------------------------------------------

CLUSTAL/aln
===========

The CLUSTAL family of multiple sequence alignment programs produce 'aln'
format.

MView option: ``-in clustal``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `CLUSTAL`    1.60, 1.70, 1.83, 2.1                  ok
 ===========  ====================================== ============


HSSP
====

The HSSP multiple sequence alignment format is produced by the MaxHom protein
sequence and structure homology alignment program.

MView option: ``-in hssp``

 ===========  ====================================== ============
 Program      Tested                                 Status
 ===========  ====================================== ============
 `HSSP`       1.0 1991                               ok
 ===========  ====================================== ============


MAF
===

The UCSC Multple Alignment Format.

MView option: ``-in maf``

MAF files contain distinct alignment blocks. By default, only the first one
will be output. You can extract a particular block, say the second one, with
``-block 2``, or all blocks with ``-block '*'``.


FASTA
=====

The classic FASTA (or Pearson) sequence file format.

MView option: ``-in fasta``


MSF
===

The Wisconsin Package GCG Multiple Sequence File format.

MView option: ``-in msf``


PIR
===

The Protein Information Resource sequence file format.

MView option: ``-in pir``


plain
=====

A simple multiple sequence alignment format.

MView option: ``-in plain``

This is composed of rows of identifier and sequence in two columns like:

.. raw:: html

  <PRE>
  identifier1   sequence1
  identifier2   sequence2
  identifier3   sequence3
  </pre>

and can contain comment lines starting with ``#``. Identifiers and sequences
must not contain any whitespace as this is used to separate the columns. The
sequences need not be aligned vertically, but they must all be the same
length. Use ``-`` and/or ``.`` characters for gaps.


Unsupported
===========

A few other formats were implemented for specific use-cases and are not
maintained:

 =============  =============== ============
 Format         MView option	Status
 =============  ===============	============
 MIPS-ALN       ``-in mips``	experimental
 MULTAS/MULTAL  ``-in multas``	experimental
 jnet -z        ``-in jnet``	experimental
 =============  ===============	============

.. END
