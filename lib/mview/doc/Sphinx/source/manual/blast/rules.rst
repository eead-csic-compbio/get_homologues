
.. _ref_blast_rules:

BLAST HSP processing rules
==========================

A BLAST run comprises, among other sections, a 'ranking' of hits with one-line
summaries including scoring information, followed by the search database
'hits' themselves. Each hit may contain one or more HSPs (ungapped BLAST) or
gapped alignments (gapped BLAST), not all of which may be reported in the
ranking.


Output modes
------------

Some control over which HSPs are to be processed by MView is provided by the
command line option ``-hsp``. There are three choices: 

 ================= ============================================================
 MView option      Description
 ================= ============================================================
 ``-hsp ranked``   Only those HSPs or alignments that contributed to the
                   ranked hit are selected. These are tiled into a single
		   alignment row per hit. This is the default behaviour.

 ``-hsp discrete`` All HSPs or alignments are selected, but each is assigned
                   its own row in the resulting stacked alignment. This is way
                   to view or extract every distinct HSP.

 ``-hsp all``      All HSPs or alignments are selected and tiled into a single
                   alignment row per hit. This is probably not what you want
                   and the option may be removed in future.
 ================= ============================================================

Precise application of these processing modes varies slightly with each BLAST
program because properties such as query and/or hit nucleotide strand
orientation or reading frame must be taken into consideration.


Algorithm for tiling HSPs
-------------------------

The ranking contains scoring information for the best alignment fragment
comprising a score (score or bits, depending on BLAST version), an estimate of
significance (E-value or p-value, depending on version) and an optional
fragment count, N.

Each database search hit comprises 1 or more fragments or HSPs, each of which
has a corresponding score, significance level and the fragment count, N. These
are grouped by significance and fragment count and are generally ordered by
decreasing score (best first).

Very often there may be just one alignment per hit and this is extracted and
output by MView. Otherwise, in ranked mode (the default), MView applies
the following algorithm to choose HSPs to tile.

**Select the ranked HSPs**

+ Read off the blast score and N from the ranking, guessing N=1 if absent.
+ Find the corresponding HSP, which should be the first one in the list of
  aligned fragments and should match the query orientation (BLASTN,
  etc). Differences in numeric rounding of the score as well as missing N
  values must be allowed for.
+ Read off the corresponding significance score for that HSP; usually it
  matches the one stated in the ranking.
+ Search the remainder of the HSPs looking for the rest of the group with the
  same N, significance and query and hit orientations.
+ Return the group of correlated HSPs.

**Select a sequentially coherent subset**

+ Take the initially found HSP with the highest score as a seed of of the
  'current set' and record its query and hit coordinate midpoints.
+ Compute the midpoints of each remaining HSP in both query sequence and
  hit sequence coordinates.
+ Sort the HSP set by midpoint in query order.
+ Search the list linearly for successive HSPs that lie upstream of the
  current set in both query and hit sequence ordering, appending these to the
  growing current set; discard HSPs that violate either query or hit sequence
  ordering. Update the current set query and hit midpoints.
+ Reverse the HSP list.
+ Search the list linearly for successive HSPs that lie downstream in both
  query and hit sequence ordering, appending these to the growing set; discard
  HSPs that violate either query or hit sequence ordering. Update the current
  set query and hit midpoints.

**Strip query gaps**

+ Excise columns in each sequence of a fragment pair to remove query gaps.
  These are generally marked for display by lowercasing the bounding symbols.

**Tile**

+ Sort the new HSP set best first by decreasing alignment score and increasing
  length.
+ Paint successive HSPs onto a pair of growing sequence scaffolds (one for the
  query, common to all alignment rows; one for the hit), using the query
  sequence coordinates as guides. All rows therefore potentially contribute to
  filling in the query sequence.
+ Resolve clashes of symbols at any position by applying a write-once policy,
  i.e., assume the best HSPs are those that are painted earlier.

Clashes lead to small discrepancies in the affected sequences and occur when
overlapping fragments have:

- slightly different alignments and/or gap positions (may affect all BLAST
  programs if the fragment count N > 1);
- different reading frames (affects BLASTX, TBLASTN, TBLASTX).

The ``-hsp discrete`` mode can be used to output each HSP in its own alignment
row, and will report all HSPs, including the alternative ones with other
significance scores.

Lastly, if you need a real sequence alignment, use the information provided by
MView to choose sequence identifiers then use a proper multiple alignment tool
such as CLUSTAL on the original sequences.


More details about each BLAST program
-------------------------------------

A more detailed description of the selection rules is given for each BLAST
family and program, as follows:


.. toctree::

   blast1
   blast2


.. END
