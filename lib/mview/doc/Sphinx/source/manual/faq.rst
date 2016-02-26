Frequently asked questions
^^^^^^^^^^^^^^^^^^^^^^^^^^


Can MView extract BLAST (or other) output into a FASTA/PIR/MSF/CLUSTAL file?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Yes - MView extracts and converts results into FASTA, PIR, MSF and CLUSTAL/aln
formats. See :ref:`ref_output_formats`.


Can MView extract or view the separate BLAST alignment or HSPs?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Yes - MView normally displays or extracts only the top scoring alignments or
HSPs for the query. You can get all alignments and/or HSPs using the ``-hsp
discrete`` option.


Does MView work with nucleotide sequences?
""""""""""""""""""""""""""""""""""""""""""

Yes - use one of the ``-moltype dna``, ``-moltype rna`` or ``-moltype na``
options to tell MView to use an appropriate colour scheme for HTML output
(currently, these three options are equivalent).


Can MView use CLUSTAL colours?
""""""""""""""""""""""""""""""

Yes - use the ``-colormap CLUSTAL`` option to change the colourmaps for
both protein and nucleotide alignments.


Can MView process data from a Web page?
"""""""""""""""""""""""""""""""""""""""

Basically, no, unless you are lucky or prepared to edit the Web page. The
MView parsers are all built to recognise the raw text output produced by the
respective programs (BLAST, FASTA, etc.) or to recognise particular flat-file
formats (MSF, PIR, etc.). When a site adds HTML markup to this to make a Web
page, arbitrary parts are changed/deleted/added polluting the text, so that
even dumping the page in text-only format still leaves traces.


Can I switch off HTML markup?
"""""""""""""""""""""""""""""

Yes - the program defaults to plain text output unless the ``-html`` option is
set to any valid value other than ``off``, or the output format is set to PIR,
MSF, FASTA, or RDB.


How can I print?
""""""""""""""""

From the web browser. To produce something that fits on typical paper sizes
one must set ``-width 60`` or similar and turn off some of the leading text,
e.g., ``-label2 -label3``.


Can I edit the alignment?
"""""""""""""""""""""""""

No - MView isn't an alignment editor. You could try to copy/paste the output
into a spreadsheet and process it there, then reload it into MView.


What do the percent identities mean?
""""""""""""""""""""""""""""""""""""

By default, percent identities reported in each alignment row are calculated
with respect to the aligned portion of the reference sequence (usually the
query or first row) thus:

.. math::

  \frac{\mathrm{number of identical residues}}
       {\mathrm{length of ungapped reference sequence over aligned region}}
  \times 100

Still, in the case of BLAST MView output, minor deviations from the
percentages reported by BLAST are due to (1) different rounding, and (2) the
way MView assembles a single pseudo-sequence (see :ref:`ref_funny_sequences`)
for a hit composed of multiple HSPs, giving an averaged percent identity. This
default behaviour above is also obtained using the option ``-pcid aligned``.

Two other calculation possibilities are available: ``-pcid reference``
normalises by the ungapped length of the query or reference sequence, and
``-pcid hit`` normalises by the ungapped length of the hit sequence.

You can change the reference sequence against which identities are calculated
using the ``-reference`` option, which requires either a row number or a
sequence identifier.


.. _ref_funny_sequences:

Why are some symbols lowercased?
""""""""""""""""""""""""""""""""

Gapped input (e.g., FASTA, BLAST2, PSI-BLAST) is subject to a further
processing step when producing the stacked alignment. The query sequence acts
as a template, but any gaps introduced into the query by a database hit are
excised to ensure a contiguous query string. The position of the excision in
the affected hit is marked by lowercasing the pair of boundary symbols.


Why is the query sequence incomplete?
"""""""""""""""""""""""""""""""""""""

The displayed query sequence is assembled from the input data, not by
reference to any external sequence database. MView will pad missing query
sequence with 'X' characters based on the numeric match ranges if this is
needed to complete an assembly. Occasionally, you may see a '?' character -
this means that a non-standard residue was seen on input.


How are overlapping BLAST HSPs processed?
"""""""""""""""""""""""""""""""""""""""""

Ungapped BLAST input is processed to produce a stack of hit sequence strings
aligned against a contiguous query sequence. The query sequence acts as a
template for each hit sequence onto which hit fragments are overlayed in the
query positions.

In outline the default method of processing of HSPs is as follows:

For BLAST (series 1), as of MView version 1.37, only the HSPs contributing to
the ranked hit contribute to this overlay process. A sorting scheme ensures
that the best of these fragments are overlayed last and are not obscured by
weaker ones, for example, BLAST hits are sorted by score and
length. Differences of ordering of fragments along query and hit naturally
result in a patchwork that may not correspond exactly to the real hit
sequences. Nevertheless, the resulting alignment stack is very informative,
and the user can always run and view a gapped search if that is preferred.

For BLAST (series 2) and PSI-BLAST, often only a single gapped alignment is
reported by blast for a given database hit. However, sometimes there are
alternative alignments and the same stacking rules apply.

Greater control over the choice of HSPs is available through the ``-hsp``
option. See :ref:`ref_blast_rules` for more details.


.. END
