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

 ========== ======= ======= ======= ====== ====== ====== ====== ====== ====== ====== =======
 Program    2.2.25+ 2.2.28+ 2.2.31+ 2.3.0+ 2.4.0+ 2.5.0+ 2.6.0+ 2.7.1+ 2.8.1+ 2.9.0+ 2.10.0+
 ========== ======= ======= ======= ====== ====== ====== ====== ====== ====== ====== =======
 `blastp`   ok      ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `blastn`   ok      ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `blastx`   .       ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `tblastn`  .       ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `tblastx`  .       ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `psiblast` .       ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `phiblast` .       .       .       .      .      .      .      .      .      ok     ok
 ========== ======= ======= ======= ====== ====== ====== ====== ====== ====== ====== =======

MView option: ``-in blast`` (with `-outfmt 7` tabular BLAST format)

 ========== ======= ======= ====== ====== ====== ====== ====== ====== ====== =======
 Program    2.2.28+ 2.2.31+ 2.3.0+ 2.4.0+ 2.5.0+ 2.6.0+ 2.7.1+ 2.8.1+ 2.9.0+ 2.10.0+
 ========== ======= ======= ====== ====== ====== ====== ====== ====== ====== =======
 `blastp`   ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `blastn`   ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `blastx`   ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `tblastn`  ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `tblastx`  ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `psiblast` ok      ok      ok     ok     ok     ok     ok     ok     ok     ok
 `phiblast` .       .       .      .      .      .      .      .      ok     ok
 ========== ======= ======= ====== ====== ====== ====== ====== ====== ====== =======

The main `-outfmt 7` fields or columns recognised by MView are:

- `std` to report standard fields: `qaccver`, `saccver`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`;
- `qseq`, `sseq` to report the query and sbjct sequences;
- `stitle` or `salltitles` to report the sbjct descriptions;
- `qframe` to report the query frame;
- `sframe` to report the sbjct frame.

When running blast with tabular output, minimal field specifiers suitable for
input to MView are:

    `blastp   -outfmt '7 std qseq sseq stitle'`

    `blastn   -outfmt '7 std qseq sseq stitle'`

    `blastx   -outfmt '7 std qseq sseq stitle qframe'`

    `tblastn  -outfmt '7 std qseq sseq stitle sframe'`

    `tblastx  -outfmt '7 std qseq sseq stitle qframe sframe'`

    `psiblast -outfmt '7 std qseq sseq stitle'`

Additional fields that are currently recognised by MView are:

- `staxid`
- `ssciname`
- `scomname`
- `sblastname`
- `sskingdom`
- `staxids`
- `sscinames`
- `scomnames`
- `sblastnames`
- `sskingdoms`

If present, these are passed through in the format conversions (`-out fasta`,
`-out pir`, `-out rdb`) and to the normal MView display mode where they appear
as extra columns at the right margin of the alignment. Having very variable
field widths they may appear badly aligned vertically. To disable them, use
the `-label8` option.

MView will silently ignore any fields other than those listed above.

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

 ========== ===== ===== =====
 Program    2.2.5 2.2.6 2.2.9
 ========== ===== ===== =====
 `blastp`   ok    ok    ok
 `blastn`   .     ok    ok
 `blastx`   .     ok    ok
 `tblastn`  .     ok    ok
 `tblastx`  .     ok    ok
 `psiblast` .     ok    ok
 `phiblast` .     ok    ok
 ========== ===== ===== =====


BLAST series 2.0
----------------

MView option: ``-in blast``

 ========== ===== ===== ===== ===== ===== ====== ======
 Program    2.0.2 2.0.4 2.0.5 2.0.6 2.0.9 2.0.10 2.0.14
 ========== ===== ===== ===== ===== ===== ====== ======
 `blastp`   .     ok    ok    .     ok    ok     .
 `blastn`   .     ok    ok    .     ok    .      ok
 `blastx`   .     .     ok    .     ok    .      .
 `tblastn`  .     .     ok    .     .     ok     .
 `tblastx`  .     .     ok    .     .     .      .
 `psiblast` ok    ok    ok    ok    .     ok     .
 `phiblast` .     .     .     .     ok    .      .
 ========== ===== ===== ===== ===== ===== ====== ======



BLAST series 1.4
----------------

MView option: ``-in blast``

 ========== ===== =====
 Program    1.4.7 1.4.9
 ========== ===== =====
 `blastp`   ok    ok
 `blastn`   .     ok
 `blastx`   .     ok
 `tblastn`  .     ok
 `tblastx`  .     ok
 ========== ===== =====


------------------------------------------------------------------------------

WashU-BLAST family
==================

Sequence database search programs.


WU-BLAST series 2.0
-------------------

MView option: ``-in blast``

 ========= ====== ====== ===
 Program   2.0a13 2.0a19 2.0
 ========= ====== ====== ===
 `blastp`  ok     ok     ok
 `blastn`  .      ok     ok
 `blastx`  .      ok     ok
 `tblastn` .      ok     ok
 `tblastx` .      ok     ok
 ========= ====== ====== ===

------------------------------------------------------------------------------

Uni. Virginia FASTA family
==========================

Sequence database search programs:


FASTA series 36
---------------

MView option: ``-in uvfasta``

 ========== ===== ====== ======= ====== ======= ======= ======= ======= =======
 Program    36.07 36.3.3 35.3.5a 36.3.6 36.3.7b 36.3.8e 36.3.8f 36.3.8g 36.3.8h
 ========== ===== ====== ======= ====== ======= ======= ======= ======= =======
 `fasta`    ok    ok     .       ok     ok      ok      ok      ok      ok
 `fastx`    .     ok     .       ok     ok      ok      ok      ok      ok
 `fasty`    .     .      .       ok     ok      ok      ok      ok      ok
 `tfastx`   .     ok     .       ok     ok      ok      ok      ok      ok
 `tfasty`   .     ok     .       ok     ok      ok      ok      ok      ok
 `ssearch`  .     ok     .       ok     ok      ok      ok      ok      ok
 `ggsearch` .     ok     .       ok     ok      ok      ok      ok      ok
 `glsearch` .     ok     .       ok     ok      ok      ok      ok      ok
 `fastm`    .     .      expt    expt   expt    expt    expt    expt    expt
 `fasts`    .     .      .       expt   expt    expt    expt    expt    expt
 `fastf`    .     .      .       expt   expt    expt    expt    expt    expt
 `tfastm`   .     .      .       .      expt    expt    expt    expt    expt
 `tfasts`   .     .      .       .      expt    expt    expt    expt    expt
 `tfastf`   .     .      .       .      expt    expt    expt    expt    expt
 ========== ===== ====== ======= ====== ======= ======= ======= ======= =======

where 'expt' indicates that MView parses and processes the format, but the
results may not be particularly useful, hence 'experimental'.


FASTA series 35
---------------

MView option: ``-in uvfasta``

 ========== =====
 Program    35.04
 ========== =====
 `fasta`    ok
 `tfastx`   ok
 `ssearch`  ok
 `ggsearch` ok
 `glsearch` ok
 ========== =====


FASTA series 34
---------------

MView option: ``-in uvfasta``

 ========== =======
 Program    34.26.5
 ========== =======
 `fasta34`  ok
 ========== =======


FASTA series 3.0 - 3.4
----------------------

MView option: ``-in uvfasta``

 ========== ====== ====== ====== ====== ====== ====== ====== ====== ====== ====== ======
 Program    3.0t76 3.0t82 3.1t07 3.2t01 3.2t05 3.2t07 3.2t05 3.2t07 3.3t01 3.3t07 3.4t23
 ========== ====== ====== ====== ====== ====== ====== ====== ====== ====== ====== ======
 `fasta`    ok     .      ok     ok     ok     ok     ok     ok     ok     ok     ok
 `fastx`    .      .      .      .      .      .      .      .      .      .      ok
 `fasty`    .      .      .      .      .      .      .      .      .      .      ok
 `tfasta`   .      .      .      .      .      .      .      .      .      .      ok
 `tfastx`   .      ok     .      .      .      .      .      .      .      .      ok
 `tfasty`   .      .      .      .      .      .      .      .      .      .      ok
 `tfastxy`  .      .      ok     .      .      .      .      .      .      .      .
 ========== ====== ====== ====== ====== ====== ====== ====== ====== ====== ====== ======


FASTA series 2
--------------

MView option: ``-in uvfasta``

 ========== ==== ======
 Program    2.0u 2.0u63
 ========== ==== ======
 `fasta`    ok   .
 `tfastx`   .    ok
 ========== ==== ======


FASTA series 1
--------------

MView option: ``-in uvfasta``

 ========== ======
 Program    1.6c24
 ========== ======
 `fasta`    ok
 ========== ======


------------------------------------------------------------------------------

CLUSTAL/aln
===========

The CLUSTAL family of multiple sequence alignment programs produce 'aln'
format.

MView option: ``-in clustal``

 ========= ==== ==== ==== ===
 Version   1.60 1.70 1.83 2.1
 ========= ==== ==== ==== ===
 `CLUSTAL` ok   ok   ok   ok
 ========= ==== ==== ==== ===


HSSP/Maxhom
===========

The HSSP multiple sequence alignment format is produced by the MaxHom protein
sequence and structure homology alignment program.

MView option: ``-in hssp``

 ======= ========
 Version 1.0 1991
 ======= ========
 `HSSP`  ok
 ======= ========


MAF
===

The UCSC Multple Alignment Format.

MView option: ``-in maf``

MAF files contain distinct alignment blocks. By default, all blocks will be
output. You can extract a particular block, say the second one, with ``-block
2``, or several with ``-block '1,2,7'``.


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
