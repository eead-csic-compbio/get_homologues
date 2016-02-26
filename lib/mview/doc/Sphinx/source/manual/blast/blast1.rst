NCBI BLAST series 1 and WU-BLAST series 2
-----------------------------------------

The rules used in MView for selecting and tiling HSPs for:

* NCBI BLAST series 1
* WU-BLAST series 2

BLASTP
^^^^^^
     
Searches a protein query sequence against a protein sequence database.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | pro | pro | pro | pro |  \+ |  \+ |  forward |   no     |   none       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation is always '\+' giving one alignment.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, HSPs corresponding to those in the
                   ranking are selected for assembly having: (i) same-sized
                   fragment set, *N*, (ii) same p-value, *P*. Row data
                   fields: (i) rank as integral *hitnum*, (ii)
		   identifier, (iii) description, (iv) ranked *score*, (v)
		   ranked *P(N)*, (vi) ranked *N*.

 ``-hsp all``      For each ranked hit, all HSPs are selected for
                   assembly. Row data fields: (i) rank as integral
		   *hitnum*, (ii) identifier, (iii) description, (iv) highest
		   *score*, (v) lowest *P*, (vi) total number of HSPs, *N*.

 ``-hsp discrete`` For each ranked hit, every HSP is assigned its own row.
                   Row data fields: (i) HSP rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   HSP *score*, (v) HSP *P*, (vi), *N=1*.
 ================= ===========================================================


BLASTN
^^^^^^

Searches a nucleotide query sequence against a nucleotide sequence database.

 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 |   Input   |   Output  |   Orient  | Assemble | Renumber | Filter       |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | qry | hit | qry | hit | qry | hit |          |          |              |
 +=====+=====+=====+=====+=====+=====+==========+==========+==============+
 | dna | dna | dna | dna |  \+ |  \+ |  forward |   no     | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | dna | dna |  \+ | \-- |  forward |   no     | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | dna | dna | \-- |  \+ |  reverse |   no     | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | dna | dna | \-- | \-- |  reverse |   no     | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation can be '\+' or '\--' giving two distinct alignments.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, HSPs corresponding to those in the
                   ranking are selected for assembly having: (i) same-sized
                   fragment set, *N*, (ii) same p-value, *P*, (iii) same hit
                   orientation. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) ranked
                   *score*, (v) ranked *P(N)*, (vi) ranked *N*, (vii) hit
		   orientation.

 ``-hsp all``      For each ranked hit, HSPs are selected for assembly
                   having: (i) same hit orientation. Row data fields: (i) rank
                   as compound *hitnum.alignum*, (ii) identifier, (iii)
                   description, (iv) highest *score*, (v) lowest *P*, (vi)
                   total number of HSPs, *N*, (vii) *hit* orientation.
                
 ``-hsp discrete`` For each ranked hit, every HSP is assigned its own row.
                   Row data fields: (i) HSP rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   HSP *score*, (v) HSP *P*, (vi) *N=1*, (vii) hit
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
 | dna | pro | pro | pro | \+F |  \+ |  forward |   yes    | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | pro | pro | pro |\--F |  \+ |  reverse |   yes    | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation can be '\+' or '\--' giving two distinct alignments by
orientation. Query reading frames of the same sign are merged unless ``-hsp
discrete`` is used.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, HSPs corresponding to those in the
                   ranking are selected for assembly having: (i) same-sized
                   fragment set, *N*, (ii) same p-value, *P*, (iii) same query
                   orientation. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) ranked
                   *score*, (v) ranked *P(N)*, (vi) ranked *N*, (vii) query
                   orientation.

 ``-hsp all``      For each ranked hit, HSPs are selected for assembly
                   having: (i) same query orientation. Row data fields: (i)
                   rank as integral *hitnum*, (ii) identifier, (iii)
                   description, (iv) highest *score*, (v) lowest *P*, (vi)
                   total number of HSPs, *N*, (vii) query orientation.

 ``-hsp discrete`` For each ranked hit, every HSP is assigned its own row.
                   Row data fields: (i) HSP rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   HSP *score*, (v) HSP *P*, (vi) *N=1*, (vii) query
		   orientation.
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
 | pro | dna | pro | pro |  \+ | \+F |  forward | no       | none         |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | pro | dna | pro | pro |  \+ |\--F |  forward | no       | none         |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation is always '\+' giving one alignment. Hit reading frames
of the same sign are merged unless ``-hsp discrete`` is used.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, HSPs corresponding to those in the
                   ranking are selected for assembly having: (i) same-sized
                   fragment set, *N*, (ii) same p-value, *P*, (iii) same hit
                   orientation. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) ranked
                   *score*, (v) ranked *P(N)*, (vi) ranked *N*, (vii) hit
		   orientation.

 ``-hsp all``      For each ranked hit, HSPs are selected for assembly
                   having: (i) same hit orientation. Row data fields: (i) rank
                   as compound *hitnum.alignum*, (ii) identifier, (iii)
                   description, (iv) highest *score*, (v) lowest *P*, (vi)
                   total number of HSPs, *N*, (vii) hit orientation.

 ``-hsp discrete`` For each ranked hit, every HSP is assigned its own row.
                   Row data fields: (i) HSP rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   HSP *score*, (v) HSP *P*, (vi) *N=1*, (vii) hit
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
 | dna | dna | pro | pro | \+F | \+F | forward  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | pro | pro | \+F |\--F | forward  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | pro | pro |\--F | \+F | reverse  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+
 | dna | dna | pro | pro |\--F |\--F | reverse  | yes      | query orient |
 +-----+-----+-----+-----+-----+-----+----------+----------+--------------+

The query orientation can be '\+' or '\--' giving two distinct alignments by
orientation. Query and hit reading frames of the same sign are merged,
respectively, unless ``-hsp discrete`` is used.

 ================= ===========================================================
 MView option      Description
 ================= ===========================================================
 ``-hsp ranked``   For each ranked hit, HSPs corresponding to those in the
                   ranking are selected for assembly having: (i) same-sized
                   fragment set, *N*, (ii) same p-value, *P*, (iii) same hit
                   orientation. Row data fields: (i) rank as integral
                   *hitnum*, (ii) identifier, (iii) description, (iv) ranked
                   *score*, (v) ranked *P(N)*, (vi) ranked *N*, (vii) hit
                   orientation.

 ``-hsp all``      For each ranked hit, HSPs are selected for assembly
                   having: (i) same hit orientation. Row data fields: (i) rank
                   as compound *hitnum.alignum*, (ii) identifier, (iii)
                   description, (iv) highest *score*, (v) lowest *P*, (vi)
                   total number of HSPs, *N*, (vii) hit orientation.

 ``-hsp discrete`` For each ranked hit, every HSP is assigned its own row.
                   Row data fields: (i) HSP rank as compound rank
                   *hitnum.alignum*, (ii) identifier, (iii) description, (iv)
		   HSP *score*, (v) HSP *P*, (vi) *N=1*, (vii) hit
		   orientation.
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
                properties, such as BLASTN query strand orientation, etc.
 ============== ==============================================================


Key
^^^

 ============== ==============================================================
 Term           Description
 ============== ==============================================================
 *hitnum*       The ordinal number (starting at 1) of a search hit as
                encountered in the ranking section of the raw BLAST output.

 *alignum*      The ordinal number (starting at 1) of a single HSP within a
                given hit (starting afresh with each hit) as encountered in the
                raw BLAST output.

 *N*            A count of HSPs which is context dependent in MView. In raw
                BLAST output it is the number of HSPs deemed by BLAST to
		belong in a coherent set of alignment fragments, but the MView
		HSP processing options (``ranked``, ``all``, ``discrete``)
		change the precise meaning.

 *P(N)*         The BLAST alignment p-value associated with a set of *N* HSPs.

 *P*            BLAST 1.4 series reports relatively large p-vales in floating
                point format rather than scientific notation. The original HSP
                p-value may be rounded to 2 decimal places in the hit ranking,
                so comparison of these values by MView when selecting HSPs
                takes this into account, which could lead to over-selection of
                HSPs on occasion.

 *score*        The BLAST alignment score.
 ============== ==============================================================


Differences between NCBI BLAST 1 and NCBI BLAST 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* In the default gapped mode, the BLAST 2 programs produce no *N* fragment
  counts. Fragment counts are switched on in BLAST when it is run in ungapped
  mode using the ``blastall... -g F`` command line option.

* Significances are reported only as Expectation (E-values), *E* rather than
  p-values, *P*.

* BLASTX: query strands don't carry frame information (early versions only).

* TBLASTN: hit strands don't carry frame information (early versions only).

* TBLASTX: query and hit strands don't carry frame information (early versions
  only).

With BLAST 2 series output, MView assumes that only one fragment (a gapped
alignment) is required in 'ranked' mode and accepts the FIRST alignment with
score and expectation matching those in the hit ranking. If you wish to see
the effect of any other HSPs use the ``-hsp all`` or ``-hsp discrete`` options
to MView. This is particularly important if you use BLAST in ungapped mode.


.. END
