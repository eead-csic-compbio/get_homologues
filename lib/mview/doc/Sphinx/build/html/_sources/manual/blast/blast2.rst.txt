NCBI BLAST series 2 and NCBI BLAST+
-----------------------------------

The rules used in MView for selecting and tiling HSPs for:

* NCBI BLAST series 2 (including PSI-BLAST)
* NCBI BLAST\+


BLASTP
^^^^^^

Searches a protein query sequence against a protein sequence database.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | pro | pro | pro | pro | \+  | \+  | forward  | no       | none         |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation is always '\+' giving one alignment.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, alignments corresponding to those in
                   the ranking are selected for assembly having: (i) same
                   expectation, *E*, (ii) same score, *bits*. Row data
                   fields: (i) rank as integral *hitnum*, (ii)
                   identifier, (iii) description, (iv) ranked *bits*, (v)
                   ranked *E*, (vi) *N=1*, (vii) hit orientation.

 ``-hsp all``      For each ranked hit, all alignments are selected for
                   assembly. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) highest
                   *bits*, (v) lowest *E*, (vi) total number of alignments
		   *N*.

 ``-hsp discrete`` For each ranked hit, every alignment is assigned its own
                   row. Row data fields: (i) alignment rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   alignment *bits*, (v) alignment *E*, (vi) *N=1*.
 ================= ===========================================================


BLASTN
^^^^^^

Searches a nucleotide query sequence against a nucleotide sequence database.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | dna | dna | dna | dna | \+  | \+  | forward  | no       | none         |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | dna | dna | \+  | \-- | forward  | no       | none         |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation can be '\+' or '\--' giving two distinct alignments.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, alignments corresponding to those in
                   the ranking are selected for assembly having: (i) same
                   expectation, *E*, (ii) same score, *bits*. Row data
                   fields: (i) rank as integral *hitnum*, (ii)
                   identifier, (iii) description, (iv) ranked *bits*, (v)
                   ranked *E*, (vi) *N=1*, (vii) hit orientation.

 ``-hsp all``      For each ranked hit, alignments are selected for assembly
                   having: (i) same hit orientation. Row data fields: (i) rank
                   as compound *hitnum.alignum*, (ii) identifier, (iii)
                   description, (iv) highest *bits*, (v) lowest *E*, (vi)
		   total number of alignments *N*, (vii) hit orientation.

 ``-hsp discrete`` For each ranked hit, every alignment is assigned its own
                   row. Row data fields: (i) alignment rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   alignment *bits*, (v) alignment *E*, (vi) *N=1*, (vii) hit
		   orientation.
 ================= ===========================================================


BLASTX
^^^^^^

Searches a nucleotide query sequence against a protein sequence database,
translating the query in all 6 reading frames.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | dna | pro | pro | pro | \+  | \+  | forward  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | pro | pro | pro | \-- | \+  | reverse  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation can be '\+' or '\--' giving two distinct alignments by
orientation. Query reading frames of the same sign are merged unless ``-hsp
discrete`` is used.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, alignments corresponding to those in
                   the ranking are selected for assembly having: (i) same
                   expectation, *E*, (ii) same score, *bits*, (iii) same query
                   orientation. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) ranked
                   *bits*, (v) ranked *E*, (vi) *N=1*, (vii) query
		   orientation.

 ``-hsp all``      For each ranked hit, alignments are selected for assembly
                   having: (i) same query orientation. Row data fields: (i)
                   rank as integral *hitnum*, (ii) identifier, (iii)
                   description, (iv) highest *bits*, (v) lowest *E*, (vi)
		   total number of alignments, *N*, (vii) query orientation.

 ``-hsp discrete`` For each ranked hit, every alignment is assigned its own
                   row. Row data fields: (i) alignment rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   alignment *bits*, (v) alignment *E*, (vi) *N=1*, (vii)
		   query orientation.
 ================= ===========================================================


TBLASTN
^^^^^^^

Searches a protein query sequence against a nucleotide sequence database,
translating the hits in all 6 reading frames.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | pro | dna | pro | pro | \+  | \+  | forward  | no       | none         |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | pro | dna | pro | pro | \+  | \-- | forward  | no       | none         |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation is always '\+' giving one alignment. Hit reading frames
of the same sign are merged unless ``-hsp discrete`` is used.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, alignments corresponding to those in
                   the ranking are selected for assembly having: (i) same
                   expectation, *E*, (ii) same score, *bits*, (iii) same hit
                   orientation. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) ranked
                   *bits*, (v) ranked *E*, (vi) *N=1*, (vii) hit orientation.

 ``-hsp all``      For each ranked hit, alignments are selected for assembly
                   having: (i) same hit orientation. Row data fields: (i) rank
                   as compound *hitnum.alignum*, (ii) identifier, (iii)
                   description, (iv) highest *bits*, (v) lowest *E*, (vi)
		   total number of alignments, *N*, (vii) hit orientation.

 ``-hsp discrete`` For each ranked hit, every alignment is assigned its own
                   row. Row data fields: (i) alignment rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   alignment *bits*, (v) alignment *E*, (vi) *N=1*, (vii) hit
		   orientation.
 ================= ===========================================================


TBLASTX
^^^^^^^

Searches a nucleotide query sequence against a nucleotide sequence database,
translating both query and hit in all 6 reading frames.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | dna | dna | pro | pro | \+  | \+  | forward  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | pro | pro | \+  | \-- | forward  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | pro | pro | \-- | \+  | reverse  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | pro | pro | \-- | \-- | reverse  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation can be '\+' or '\--' giving two distinct alignments by
orientation. Query and hit reading frames of the same sign are merged,
respectively, unless ``-hsp discrete`` is used.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, alignments corresponding to those in
                   the ranking are selected for assembly having: (i) same
                   expectation, *E*, (ii) same score, *bits*, (iii) same hit
                   orientation. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) ranked
                   *bits*, (v) ranked *E*, (vi) *N=1*, (vii) hit orientation.

 ``-hsp all``      For each ranked hit, alignments are selected for assembly
                   having: (i) same hit orientation. Row data fields: (i) rank
                   as compound *hitnum.alignum*, (ii) identifier, (iii)
                   description, (iv) highest *bits*, (v) lowest *E*, (vi)
		   total number of alignments, *N*, (vii) hit orientation.

 ``-hsp discrete`` For each ranked hit, every alignment is assigned its own
                   row. Row data fields: (i) alignment rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   alignment *bits*, (v) alignment *E*, (vi) *N=1*, (vii) hit
		   orientation.
 ================= ===========================================================


PSI-BLAST
^^^^^^^^^

Iteratively searches a protein query sequence against a protein sequence
database. PSI-BLAST output is processed as a sequence of BLASTP cycles.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | pro | pro | pro | pro | \+  | \+  | forward  | no       | cycle        |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

BLASTP and PSI-BLAST share the same HSP processing rules. The query
orientation is always '\+' giving one alignment per cycle.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, the first alignment is selected to
                   form an alignment row using: (i) same expectation,
		   *E*, (ii) same score, *bits*. Row data fields: (i) rank as
		   integral *hitnum*, (ii) identifier, (iii) description, (iv)
		   ranked *bits*, (v) ranked *E*, (vi) *N=1*.

 ``-hsp all``      For each ranked hit, all alignments are selected for
                   assembly. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) highest
                   *bits*, (v) lowest *E*, (vi) total number of alignments
		   *N*.

 ``-hsp discrete`` For each ranked hit, every alignment is assigned its own
                   row. Row data fields: (i) alignment rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   alignment *bits*, (v) alignment *E*, (vi) *N=1*.
 ================= ===========================================================


Column headings
^^^^^^^^^^^^^^^

 ============== ==============================================================
 Heading        Description
 ============== ==============================================================
 Input          The input sequence types for the query, 'qry', and the 'hit',
                sequences which may be protein 'pro' or nucleotide 'dna'.

 Output         The output sequence types produced by the search, may be
                protein 'pro' or nucleotide 'dna'.

 Orient         The orientation of the query and hit sequences given as
                forward '\+' or reverse '\--' and sometimes with a reading
		frame number indicated by '\+F' or '\--F'.

 Assemble       The sequence fragment assembly algorithm:

                *forward*: Given a list of fragments each numbered from low
                to high along the query strand (i.e., conventional
                orientation): stack these at these positions along the query
                template.

                *reverse*: Given a list of fragments each numbered from
                high to low along the query strand (i.e., individually
                reversed): reverse each fragment, stack onto the template,
                reverse the template, which recovers the original numbering
                high to low and preserves intra- and inter-fragment order.

 Renumber       Whether the sequence numbering is modified by MView or not.
                Given a :math:`(start, stop)` pair in nucleotide units, the
                corresponding nearest values in amino acid units are given by:
                :math:`(start\prime, stop\prime) = (int((start+2)/3),
		int(stop/3))` This renumbering is applied to each fragment's
		positions along its parent query strand to determine the query
		template positions in amino acid units.

 Filter         Separate MView output must be produced for varous BLAST hit
                properties, such as PSI-BLAST search cycle, BLASTN query
		strand orientation, etc.
 ============== ==============================================================


Key
^^^

 ============== ==============================================================
 Term           Description
 ============== ==============================================================
 *hitnum*       The ordinal number (starting at 1) of a search hit as
                encountered in the ranking section of the raw BLAST output.

 *alignum*      The ordinal number (starting at 1) of a single HSP within a
                given hit (starting afresh with each hit) as encountered in
		the raw BLAST output.

 *N*            A count of HSPs which is context dependent in MView. In raw
                BLAST output it is the number of HSPs deemed by BLAST to
		belong in a coherent set of alignment fragments, but the MView
		HSP processing options (``ranked``, ``all``, ``discrete``)
		change the precise meaning. BLAST series 2 onwards produces
		gapped alignments by default rather than separate HSPs so
		*N=1*.

 *E*            Expectation or e-value: the expected number of hits by chance.

 *bits*         NCBI BLAST 2.0 hit ranking rounds HSP scores which are
                floating point to the nearset integer, so MView takes this
		into account when selecting HSPs.
 ============== ==============================================================

.. END
