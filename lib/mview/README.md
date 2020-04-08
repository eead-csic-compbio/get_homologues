### Introduction

**MView** is a command line utility that extracts and reformats the results of
a sequence database search or a multiple alignment, optionally adding HTML
markup for web page layout. It can also be used as a filter to extract and
convert searches or alignments to common formats.

Inputs:

- Sequence database search: BLAST, FASTA suites.
- Multiple sequence alignment: CLUSTAL, HSSP, MSF, FASTA, PIR, MAF.

Outputs:

- HTML, FASTA, CLUSTAL, MSF, PIR, RDB (tab-separated).

The tool is used in molecular biology and biomedical research for data
analyses and as a component in various bioinformatics web services. Research
papers citing MView are indexed on [Google
Scholar](https://scholar.google.com/citations?user=4ughzM0AAAAJ&hl=en "MView
citations").


### Manual

Full documentation can be found in the
[manual](https://desmid.github.io/mview/ "MView manual"), a copy of which is
bundled with the code.


### Requirements

MView is implemented in Perl, version 5 as a self-contained command line
program that should run cross-platform.

[Perl](https://www.perl.org/) is generally installed on Linux and UNIX
systems. MView is known to work on Windows with
[Strawberry Perl](http://strawberryperl.com/).


### Download

- The [current release](
  https://sourceforge.net/projects/bio-mview/files/bio-mview/mview-1.67/
  "MView current release on SourceForge") can be downloaded from SourceForge
  as a gzip or bzip2 compressed tar archive.
- Older [releases](
  https://sourceforge.net/projects/bio-mview/files/bio-mview/ "MView older
  releases on SourceForge") and historical [download statistics](
  https://sourceforge.net/projects/bio-mview/files/stats/timeline?dates=2005-01-01+to+2025-01-01
  "MView download statistics") can also be found on SourceForge.
- The [latest code](https://github.com/desmid/mview "MView source") can be
  downloaded direct from GitHub by clicking the green "Clone or download"
  button and following the instructions to either clone the git repository or
  download a ZIP archive.


### Installation

There are several ways to install MView:

- [Linux, Apple, UNIX](README_unix.md)
- [Windows](README_windows.md)
- [Perl module](README_perl.md) (advanced)


### Testing

Each release of MView is regression tested against hundreds of sample data
inputs for all the sequence database search and alignment formats and versions
thereof that are supported, together with known edge cases. This is well over
0.5GB of material, so it's not currently available externally.


### Found a bug?

Please open an issue on the MView [issue tracker](https://github.com/desmid/mview/issues "issue tracker") or send an email to `biomview@gmail.com`.

If MView isn't able to parse your input file or produces a warning message, it
would be very helpful if you can include/attach the data file in your message
so that I can (1) quickly reproduce the error, and (2) add the example to the
test suite.


### Citation

If you use MView in your work, please cite:

> Brown, N.P., Leroy C., Sander C. (1998). MView: A Web compatible database
> search or multiple alignment viewer. *Bioinformatics*. **14** (4):380-381.
> [[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/9632837 "PubMed link")]


### Copyright and licence

MView is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.


### Acknowledgements

People who contributed early code or suggestions include C. Leroy and other
members of the former Sander group at EBI. Useful suggestions relating to the
EBI sequence database search services have come from R. Lopez, W. Li and
H. McWilliam at EBI. Thanks to the many other people who have suggested new
features and reported bugs. Finally, thank you to everyone who has cited MView
in their publications.
