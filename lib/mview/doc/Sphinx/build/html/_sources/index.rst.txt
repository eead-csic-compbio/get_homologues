.. raw:: html

  <LINK rel="stylesheet" type="text/css" href="_static/MView.css">

.. toctree::
   :hidden:

   contents
   install/unix
   install/windows
   install/perl


MView
=====

**MView** is a command line utility that extracts and reformats the results of
a sequence database search or a multiple alignment, optionally adding HTML
markup for web page layout. It can also be used as a filter to extract and
convert searches or alignments to common formats.

.. mview -in fasta -html head -css on -ruler on -srs on -colormap clustal
..   -coloring consensus -threshold 90 -consensus on -con_threshold 90,80
..   -con_ignore class -con_coloring any data.dat

.. raw:: html


  <PRE>
                                  1 [        .         .         .         .         :         .         .         ] 80
  1 <A HREF="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=P00533">sp|P00533.2|EGFR_HUMAN</A>  100.0%  <SPAN CLASS=S38>F</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">KI</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">V</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">S</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S38>AF</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">T</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">YK</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">L</SPAN><SPAN CLASS=S38>W</SPAN><SPAN style="color:#666666">IPEGEK---------VKIP</SPAN><SPAN CLASS=S38>VAI</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S37>R</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S38>A</SPAN><SPAN style="color:#666666">TSPK-ANK</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>I</SPAN><SPAN style="color:#666666">L</SPAN><SPAN CLASS=S42>DE</SPAN><SPAN CLASS=S38>A</SPAN><SPAN style="color:#666666">Y</SPAN><SPAN CLASS=S38>VM</SPAN><SPAN style="color:#666666">A</SPAN><SPAN CLASS=S39>S</SPAN><SPAN CLASS=S38>V</SPAN><SPAN CLASS=S42>D</SPAN><SPAN CLASS=S39>N</SPAN><SPAN style="color:#666666">P</SPAN><SPAN CLASS=S40>H</SPAN><SPAN CLASS=S38>V</SPAN><SPAN CLASS=S41>C</SPAN><SPAN CLASS=S37>R</SPAN><SPAN CLASS=S38>LL</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S38>I</SPAN><SPAN CLASS=S41>C</SPAN>
  2 <A HREF="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=Q9I7F7">sp|Q9I7F7.3|PR2_DROME</A>    35.7%  <SPAN CLASS=S38>I</SPAN><SPAN CLASS=S39>S</SPAN><SPAN style="color:#666666">VN</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">Q</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">T</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>F</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">I</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">QQ</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">V</SPAN><SPAN CLASS=S38>W</SPAN><SPAN style="color:#666666">SNGNE-----------RIQ</SPAN><SPAN CLASS=S38>VAI</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S41>C</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S41>C</SPAN><SPAN style="color:#666666">R</SPAN><SPAN CLASS=S42>E</SPAN><SPAN style="color:#666666">RMQS-NPM</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>F</SPAN><SPAN style="color:#666666">L</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>A</SPAN><SPAN style="color:#666666">A</SPAN><SPAN CLASS=S38>IM</SPAN><SPAN style="color:#666666">H</SPAN><SPAN CLASS=S39>S</SPAN><SPAN CLASS=S38>I</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S40>H</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S39>N</SPAN><SPAN CLASS=S38>IV</SPAN><SPAN CLASS=S37>R</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S40>Y</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S38>VV</SPAN>
  3 <A HREF="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=Q08881">sp|Q08881.1|ITK_HUMAN</A>    32.9%  <SPAN CLASS=S38>L</SPAN><SPAN CLASS=S39>T</SPAN><SPAN style="color:#666666">FV</SPAN><SPAN CLASS=S39>Q</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S38>I</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">S</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S39>Q</SPAN><SPAN CLASS=S38>F</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">L</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">HL</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">Y</SPAN><SPAN CLASS=S38>W</SPAN><SPAN style="color:#666666">LN--------------KDK</SPAN><SPAN CLASS=S38>VAI</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S39>T</SPAN><SPAN CLASS=S38>I</SPAN><SPAN CLASS=S37>R</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">AMS---EE</SPAN><SPAN CLASS=S42>D</SPAN><SPAN CLASS=S38>F</SPAN><SPAN style="color:#666666">I</SPAN><SPAN CLASS=S42>EE</SPAN><SPAN CLASS=S38>A</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S38>VM</SPAN><SPAN style="color:#666666">M</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S39>S</SPAN><SPAN CLASS=S40>H</SPAN><SPAN style="color:#666666">P</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S38>LV</SPAN><SPAN CLASS=S39>Q</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S40>Y</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S38>V</SPAN><SPAN CLASS=S41>C</SPAN>
  4 <A HREF="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=Q13308">sp|Q13308.2|PTK7_HUMAN</A>   21.2%  <SPAN CLASS=S38>I</SPAN><SPAN CLASS=S37>R</SPAN><SPAN style="color:#666666">EV</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">Q</SPAN><SPAN CLASS=S38>I</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">V</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S39>Q</SPAN><SPAN CLASS=S38>F</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">A</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">VL</SPAN><SPAN CLASS=S38>A</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S38>M</SPAN><SPAN style="color:#666666">TGLS-XLPKGSMNADGVAL</SPAN><SPAN CLASS=S38>VAV</SPAN><SPAN CLASS=S37>KK</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">P</SPAN><SPAN CLASS=S42>D</SPAN><SPAN style="color:#666666">VSD-EVLQ</SPAN><SPAN CLASS=S39>S</SPAN><SPAN CLASS=S38>F</SPAN><SPAN style="color:#666666">D</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>I</SPAN><SPAN style="color:#666666">K</SPAN><SPAN CLASS=S38>FM</SPAN><SPAN style="color:#666666">S</SPAN><SPAN CLASS=S39>Q</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S39>Q</SPAN><SPAN CLASS=S40>H</SPAN><SPAN style="color:#666666">D</SPAN><SPAN CLASS=S39>S</SPAN><SPAN CLASS=S38>IV</SPAN><SPAN CLASS=S39>Q</SPAN><SPAN CLASS=S38>LLAI</SPAN><SPAN CLASS=S41>C</SPAN>
  5 <A HREF="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?val=P34265">sp|P34265.4|KIN31_CAEEL</A>  31.5%  <SPAN CLASS=S38>V</SPAN><SPAN CLASS=S42>E</SPAN><SPAN style="color:#666666">LT</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">K</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S38>AF</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">E</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">WK</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">K</SPAN><SPAN CLASS=S38>L</SPAN><SPAN style="color:#666666">LKILDA-------NHQPVL</SPAN><SPAN CLASS=S38>VAV</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S39>T</SPAN><SPAN CLASS=S38>A</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">L</SPAN><SPAN CLASS=S42>E</SPAN><SPAN style="color:#666666">SMTKEQIK</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>I</SPAN><SPAN style="color:#666666">M</SPAN><SPAN CLASS=S37>R</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>A</SPAN><SPAN style="color:#666666">R</SPAN><SPAN CLASS=S38>LM</SPAN><SPAN style="color:#666666">R</SPAN><SPAN CLASS=S39>N</SPAN><SPAN CLASS=S38>L</SPAN><SPAN CLASS=S42>D</SPAN><SPAN CLASS=S40>H</SPAN><SPAN style="color:#666666">I</SPAN><SPAN CLASS=S39>N</SPAN><SPAN CLASS=S38>VV</SPAN><SPAN CLASS=S37>K</SPAN><SPAN CLASS=S38>FF</SPAN><SPAN CLASS=S44>G</SPAN><SPAN CLASS=S38>VA</SPAN>
    consensus/90%                   <SPAN style="color:#666666">.......</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S38>F</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">........................</SPAN><SPAN CLASS=S38>VA</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">.................</SPAN><SPAN CLASS=S42>E</SPAN><SPAN style="color:#666666">...</SPAN><SPAN CLASS=S38>M</SPAN><SPAN style="color:#666666">...............</SPAN>
    consensus/80%                   <SPAN style="color:#666666">....</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">..</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S38>F</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">..</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">.....................</SPAN><SPAN CLASS=S38>VA</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S37>K</SPAN><SPAN style="color:#666666">.................</SPAN><SPAN CLASS=S42>E</SPAN><SPAN CLASS=S38>A</SPAN><SPAN style="color:#666666">..</SPAN><SPAN CLASS=S38>M</SPAN><SPAN style="color:#666666">....</SPAN><SPAN CLASS=S40>H</SPAN><SPAN style="color:#666666">...</SPAN><SPAN CLASS=S38>V</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S38>L</SPAN><SPAN style="color:#666666">.</SPAN><SPAN CLASS=S44>G</SPAN><SPAN style="color:#666666">..</SPAN>
  </PRE>

Inputs:

- Sequence database search: BLAST, FASTA suites.
- Multiple sequence alignment: CLUSTAL, HSSP, MSF, FASTA, PIR, MAF.

Outputs:

- HTML, FASTA, CLUSTAL, MSF, PIR, RDB (tab-separated).

The tool is used in molecular biology and biomedical research for data
analyses and as a component in various bioinformatics web services. Research
papers citing MView are indexed on `Google Scholar`_.

.. _Google Scholar:  https://scholar.google.com/citations?user=4ughzM0AAAAJ&hl=en


Documentation
^^^^^^^^^^^^^

Short contents:

* `Contents <contents.html>`_ (long)
* `Manual <manual/manual.html>`_
* `Input formats <manual/formats_in.html>`_
* `Output formats <manual/formats_out.html>`_
* `FAQ <manual/faq.html>`_


Requirements
^^^^^^^^^^^^

MView is implemented in Perl, version 5 as a self-contained command line
program that should run cross-platform.

`Perl <https://www.perl.org/>`_ is generally installed on Linux and UNIX
systems. MView is known to work on Windows with `Strawberry Perl
<http://strawberryperl.com/>`_.


Download
^^^^^^^^

The `current release`_ can be downloaded from SourceForge as a gzip or bzip2
compressed tar archive.

Older `releases`_ and historical `download statistics`_ can also be found on
SourceForge.

The `latest code`_ can be downloaded direct from GitHub by clicking the green
"Clone or download" button and following the instructions to either clone the
git repository or download a ZIP archive.

.. _current release:      https://sourceforge.net/projects/bio-mview/files/bio-mview/mview-1.67/
.. _releases:             https://sourceforge.net/projects/bio-mview/files/bio-mview/
.. _download statistics:  https://sourceforge.net/projects/bio-mview/files/stats/timeline?dates=2005-01-01+to+2025-01-01
.. _latest code:          https://github.com/desmid/mview


Installation
^^^^^^^^^^^^

There are several ways to install MView:

* `Linux, Apple, UNIX <install/unix.html>`_
* `Windows <install/windows.html>`_
* `Perl module <install/perl.html>`_ (advanced)


Found a bug?
^^^^^^^^^^^^

Please open an issue on the MView `issue tracker`_ or send an email to
*biomview@gmail.com*.

.. _issue tracker: https://github.com/desmid/mview/issues

If MView isn't able to parse your input file or produces a warning message, it
would be very helpful if you can include/attach the data file in your email so
that I can (1) quickly reproduce the error, and (2) add the example to the
test suite.


Citation
^^^^^^^^

If you use MView in your work, please cite:

    Brown, N.P., Leroy C., Sander C. (1998). MView: A Web compatible database
    search or multiple alignment viewer. *Bioinformatics*. **14** (4):380-381.
    [`PubMed <http://www.ncbi.nlm.nih.gov/pubmed/9632837?dopt=Abstract>`_]


Copyright and licence
^^^^^^^^^^^^^^^^^^^^^

MView is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either `version 2`_ of the License, or (at your option) any later
version.

.. _version 2: etc/Licence.html


Acknowledgements
^^^^^^^^^^^^^^^^

People who contributed early code or suggestions include C. Leroy and other
members of the former Sander group at EBI. Useful suggestions relating to the
EBI sequence database search services have come from R. Lopez, W. Li and
H. McWilliam at EBI. Thanks to the many other people who have suggested new
features and reported bugs. Finally, thank you to everyone who has cited MView
in their publications.

.. END
