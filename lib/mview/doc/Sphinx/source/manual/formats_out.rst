.. _ref_output_formats:

==============
Output formats
==============

Instead of the pretty display output, MView can write several useful multiple
alignment formats.

Any row or column filtering options will be respected, so you can also extract
subsets of sequences (see :ref:`ref_filtering_rows`), or particular alignment
blocks or strand oriented data from blast or similar (see the various sequence
database search input formats).


plain
=====

A simple multiple sequence alignment format.

MView option: ``-out plain``

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


FASTA
=====

The classic FASTA (or Pearson) sequence file format.

MView option: ``-out fasta`` or ``-out pearson``.


CLUSTAL/aln
===========

The CLUSTAL multiple sequence alignment format (version 2.1, as produced by
clustalx-2.1).

MView option: ``-out clustal`` or ``-out aln``.

Output includes end of row numbers and a Clustal-style ``*:.`` amino acid
conservation line: ``*`` for full column identity, and ``:`` or ``.`` for
strong and weak amino acid grouping, respectively, as defined in CLUSTAL.

For DNA or RNA sequences, if the molecule type was set to nucleic acid with
``-moltype na`` or ``dna`` or ``rna``, then the clustal conservation line will
show only the column identities.


PIR
===

The Protein Information Resource sequence file format.

MView option: ``-out pir``

This format includes a molecule type annotation. If the ``-moltype`` option is
set (one of ``aa``, ``dna``, ``rna``, ``na``), MView will try to output a
reasonable value (P1, DL, RL, XX; for protein, linear DNA, linear RNA, or
unknown, respectively). The default is protein.


MSF
===

The Wisconsin Package GCG Multiple Sequence File format.

MView option: ``-out msf``

This format includes a molecule type annotation. If the ``-moltype`` option is
set (one of ``aa``, ``dna``, ``rna``, ``na``), MView will output the right
value (P or N; for protein or nucleic acid). The default is protein.


RDB
===

MView option ``-out rdb``

A tab-separated format including scoring information. RDB is a CSV-style
format with header lines giving column name and type information (see `rdb`_).

.. _rdb: http://www.drdobbs.com/rdb-a-unix-command-line-database/199101326

.. END
