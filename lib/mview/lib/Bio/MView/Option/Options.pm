# Copyright (C) 2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::MView::Option::Options;

use Bio::MView::Option::Arguments;
use Bio::MView::Option::Types;
use vars qw($Header $Groups $Options);

###########################################################################
$Header = "usage: <PROG> [options] [file...]

Option names and parameter values can generally be abbreviated. Alternative
parameter values are listed in braces {}, followed by the default value in
square brackets [].

Some options take multiple arguments which must be supplied as a comma
separated list, like '1,8,9,10'. Subranges are allowed, so you could also
write that as '1,8:10' or even '1,8..10'. Any argument must be quoted if it
contains whitespace or a wildcard that might be expanded by the shell.

Option processing can be terminated using '--'.

";

###########################################################################
$Options = [

    ##################################################
    # hidden args
    ##################################################

    {
        'group'   => "HIDDEN",
        'options' => [
            {
                'option'  => "verbose",
                'type'    => "u_int",
                'default' => 0,
                'order'   => 0,
            },
        ],
    },

    ##################################################
    # i/o formats
    ##################################################

    {
        'group'   => 'FORMATS',
        'header'  => "Input/output formats:",
        'options' => [
            {
                'option'  => "in",
                'usage'   => "Input",
                'type'    => "format::in",
                'param'   => "infmt",
            },

            {
                'option'  => "out",
                'usage'   => "Output",
                'type'    => "format::out",
                'param'   => "outfmt",
            },
        ],
    },

    ##################################################
    # formatting: content
    ##################################################

    {
        'group'   => "CONTENT",
        'header'  => "Main formatting options:",
        'options' => [
            {
                'option'  => "ruler",
                'usage'   => "Show ruler",
                'type'    => "binary",
                'default' => "on",
            },

            {
                'option'  => "alignment",
                'usage'   => "Show alignment",
                'type'    => "binary",
                'default' => "on",
            },

            {
                'option'  => "conservation",
                'usage'   => "Show clustal conservation line",
                'type'    => "binary",
                'default' => "off",
            },

            {
                'option'  => "consensus",
                'usage'   => "Show consensus",
                'type'    => "binary",
                'default' => "off",
            },

            {
                'option'  => "width",
                'usage'   => "Paginate alignment in blocks of width",
                'type'    => "content::width",
            },
        ],
    },

    ##################################################
    # formatting: identity filtering
    ##################################################

    {
        'group'   => "IDENTITY",
        'header'  => "Percent identity calculations and filters:",
        'options' => [
            {
                'option'  => "pcid",
                'usage'   => "Compute percent identities with respect to",
                'type'    => "identity::pcid",
            },

            {
                'option'  => "reference",
                'usage'   => "Use row N or row identifier as %identity reference",
                'type'    => "string",
                'default' => "query",
                'param'   => "ref_id",
            },

            {
                'option'  => "minident",
                'usage'   => "Only report sequences with percent identity >= N compared to reference",
                'type'    => "percentage",
                'default' => 0
            },

            {
                'option'  => "maxident",
                'usage'   => "Only report sequences with percent identity <= N compared to reference",
                'type'    => "percentage",
                'default' => 100,
            },

            {
                'option'  => "sort",
                'usage'   => "Resort output by coverage or percent identity",
                'type'    => "sort::mode",
                'default' => "none",
            },
        ],
    },

    ##################################################
    # formatting: general filtering
    ##################################################

    {
        'group'   => "FILTER",
        'header'  => "General row/column filters:",
        'options' => [
            {
                'option'  => "top",
                'usage'   => "Report top N hits",
                'type'    => "filter::top",
                'param'   => "topn",
            },

            {
                'option'  => "show",
                'usage'   => "Keep rows 1..N or identifiers",
                'type'    => "string_list",
                'param'   => "keeplist",
            },

            {
                'option'  => "hide",
                'usage'   => "Hide rows 1..N or identifiers",
                'type'    => "string_list",
                'param'   => "skiplist",
            },

            {
                'option'  => "nops",
                'usage'   => "Exclude rows 1..N or identifiers from calculations",
                'type'    => "string_list",
                'param'   => "nopslist",
            },

            {
                'option'  => "range",
                'usage'   => "Display column range M:N as numbered by ruler",
                'type'    => "integer_range",
                'default' => "all",
            },
        ],
    },

    ##################################################
    # molecule type
    ##################################################

    {
        'group'   => "MOLTYPE",
        'header'  => "Molecule type:",
        'options' => [
            {
                'option'  => "moltype",
                'usage'   => "Affects coloring and format converions",
                'type'    => "moltype::type",
                'action'  => sub {
                    my ($self, $pn, $pv, $e) = @_;
                    #reset default coloring schemes for subsequent options
                    $self->update_option('ALIGNMENT::colormap',
                                         get_default_sequence_colormap($pv));
                    $self->update_option('CONSENSUS::con_colormap',
                                         get_default_consensus_colormap($pv));
                    $self->update_option('ALIGNMENT::groupmap',
                                         get_default_groupmap($pv));
                    $self->update_option('CONSENSUS::con_groupmap',
                                         get_default_groupmap($pv));  #same
                },
            },
        ],
    },

    ##################################################
    # colouring: alignment
    ##################################################

    {
        'group'   => "ALIGNMENT",
        'header'  => "Alignment coloring:",
        'options' => [
            {
                'option'  => "coloring",
                'usage'   => "Basic style of coloring",
                'type'    => "alignment::coloring",
                'param'   => "aln_coloring",
            },

            {
                'option'  => "colormap",
                'usage'   => "Name of colormap to use",
                'type'    => "alignment::colormap",
                'param'   => "aln_colormap",
            },

            {
                'option'  => "groupmap",
                'usage'   => "Name of groupmap to use if coloring by consensus",
                'type'    => "alignment::groupmap",
                'param'   => "aln_groupmap",
            },

            {
                'option'  => "threshold",
                'usage'   => "Threshold percentage for consensus coloring",
                'type'    => "percentage_ge_50",
                'default' => 70,
                'param'   => "aln_threshold",
            },

            {
                'option'  => "ignore",
                'usage'   => "Ignore singleton or class groups",
                'type'    => "alignment::ignore",
                'param'   => "aln_ignore",
            },
        ],
    },

    ##################################################
    # colouring: consensus
    ##################################################
    {
        'group'   => "CONSENSUS",
        'header'  => "Consensus coloring:",
        'options' => [
            {
                'option'  => "con_coloring",
                'usage'   => "Basic style of coloring",
                'type'    => "consensus::coloring",
            },

            {
                'option'  => "con_colormap",
                'usage'   => "Name of colormap to use",
                'type'    => "consensus::colormap",
            },

            {
                'option'  => "con_groupmap",
                'usage'   => "Name of groupmap to use if coloring by consensus",
                'type'    => "consensus::groupmap",
            },

            {
                'option'  => "con_threshold",
                'usage'   => "Consensus line thresholds",
                'type'    => "consensus::threshold",
            },

            {
                'option'  => "con_ignore",
                'usage'   => "Ignore singleton or class groups",
                'type'    => "consensus::ignore",
            },

            {
                'option'  => "con_gaps",
                'usage'   => "Count gaps during consensus computations if set to 'on'",
                'type'    => "binary",
                'default' => "on",
            },
        ],
    },

    ##################################################
    # colouring: patterns / motifs
    ##################################################

    {
        'group'   => "PATTERNS",
        'header'  => "Motif colouring:",
        'options' => [
            {
                'option'  => "find",
                'usage'   => "Find and highlight exact string or simple regular expression or ':' delimited set of patterns",
                'type'    => "string",
                'label'   => "pattern",
            },
        ],
    },

    ##################################################
    # formatting: miscellaneous / labels
    ##################################################

    {
        'group'   => "MISC_FORMATTING",
        'header'  => "Miscellaneous formatting:",
        'options' => [
            {
                'option'  => "label0",
                'usage'   => "Switch off label {0= row number}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label1",
                'usage'   => "Switch off label {1= identifier}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label2",
                'usage'   => "Switch off label {2= description}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label3",
                'usage'   => "Switch off label {3= scores}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label4",
                'usage'   => "Switch off label {4= percent coverage}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label5",
                'usage'   => "Switch off label {5= percent identity}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label6",
                'usage'   => "Switch off label {6= first sequence positions: query}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label7",
                'usage'   => "Switch off label {7= second sequence positions: hit}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "label8",
                'usage'   => "Switch off label {8= trailing fields}",
                'type'    => "flag::invert",
            },

            {
                'option'  => "gap",
                'usage'   => "Use this gap character",
                'type'    => "char",
                'default' => '-',
            },

            {
                'option'  => "sequences",
                'usage'   => "Output sequences",
                'type'    => "binary",
                'default' => "on",
            },

            {
                'option'  => "register",
                'usage'   => "Output multi-pass alignments with columns in register",
                'type'    => "binary",
                'default' => "on",
            },
        ],
    },

    ##################################################
    # HTML: page-level markup
    ##################################################

    {
        'group'   => "HTML",
        'header'  => "HTML markup:",
        'options' => [
            {
                'option'  => "html",
                'usage'   => "Controls amount of HTML markup",
                'type'    => "html::html_mode",
            },

            {
                'option'  => "bold",
                'usage'   => "Use bold emphasis for coloring sequence symbols",
                'type'    => "flag",
            },

            {
                'option'  => "css",
                'usage'   => "Use Cascading Style Sheets",
                'type'    => "html::css_mode",
                'param'   => "css1",
            },

            {
                'option'  => "title",
                'usage'   => "Page title string",
                'type'    => "string",
            },

            {
                'option'  => "pagecolor",
                'usage'   => "Page backgound color",
                'type'    => "color",
                'default' => "white",
            },

            {
                'option'  => "textcolor",
                'usage'   => "Page text color",
                'type'    => "color",
                'default' => "black",
            },

            {
                'option'  => "alncolor",
                'usage'   => "Alignment background color",
                'type'    => "color",
                'default' => "white",
            },

            {
                'option'  => "labcolor",
                'usage'   => "Alignment label color",
                'type'    => "color",
                'default' => "black",
            },

            {
                'option'  => "symcolor",
                'usage'   => "Alignment symbol default color",
                'type'    => "color",
                'default' => "#666666",
            },

            {
                'option'  => "gapcolor",
                'usage'   => "Alignment gap color",
                'type'    => "color",
                'default' => "#666666",
            },
        ],
    },

    ##################################################
    # links: SRS
    ##################################################

    {
        'group'   => "SRS_LINKS",
        'header'  => "Database links:",
        'options' => [
            {
                'option'  => "srs",
                'usage'   => "Try to use sequence database links",
                'type'    => "binary",
                'default' => "off",
            },

            {
                'option'  => "linkcolor",
                'usage'   => "Link color",
                'type'    => "color",
                'default' => "blue",
            },

            {
                'option'  => "alinkcolor",
                'usage'   => "Active link color",
                'type'    => "color",
                'default' => "red",
            },

            {
                'option'  => "vlinkcolor",
                'usage'   => "Visited link color",
                'type'    => "color",
                'default' => "purple",
            },
        ],
    },

    ##################################################
    # formats: BLAST 1
    ##################################################

    {
        'group'   => "BLAST1",
        'header'  => "NCBI BLAST (series 1), WashU-BLAST:",
        'options' => [
            {
                'option'  => "hsp",
                'type'    => "blast::hsp_mode",
            },

            {
                'option'  => "maxpval",
                'type'    => "blast::maxpval",
            },

            {
                'option'  => "minscore",
                'type'    => "blast::minscore",
            },

            {
                'option'  => "strand",
                'type'    => "search::strand",
            },

            {
                'option'  => "keepinserts",
                'type'    => "blast::keepinserts",
            },
        ],
    },

    ##################################################
    # formats: BLAST 2
    ##################################################

    {
        'group'   => "BLAST2",
        'header'  => "NCBI BLAST (series 2), BLAST+:",
        'options' => [
            {
                'option'  => "hsp",
                'type'    => "blast::hsp_mode",
            },

            {
                'option'  => "maxeval",
                'type'    => "blast::maxeval",
            },

            {
                'option'  => "minbits",
                'type'    => "blast::minbits",
            },

            {
                'option'  => "strand",
                'type'    => "search::strand",
            },

            {
                'option'  => "keepinserts",
                'type'    => "blast::keepinserts",
            },
        ],
    },

    ##################################################
    # formats: PSI-BLAST
    ##################################################

    {
        'group'   => "PSIBLAST",
        'header'  => "NCBI PSI-BLAST:",
        'options' => [
            {
                'option'  => "hsp",
                'type'    => "blast::hsp_mode",
            },

            {
                'option'  => "maxeval",
                'type'    => "blast::maxeval",
            },

            {
                'option'  => "minbits",
                'type'    => "blast::minbits",
            },

            {
                'option'  => "cycle",
                'type'    => "psiblast::cycle",
            },

            {
                'option'  => "keepinserts",
                'type'    => "blast::keepinserts",
            },
        ],
    },

    ##################################################
    # formats: FASTA (U. of Virginia)
    ##################################################

    {
        'group'   => "FASTA",
        'header'  => "FASTA (U. of Virginia):",
        'options' => [
            {
                'option'  => "minopt",
                'type'    => "fasta::minopt",
            },

            {
                'option'  => "strand",
                'type'    => "search::strand",
            },
        ],
    },

    ##################################################
    # formats: HSSP / Maxhom
    ##################################################

    {
        'group'   => "HSSP",
        'header'  => "HSSP/Maxhom:",
        'options' => [
            {
                'option'  => "chain",
                'type'    => "hssp::chain",
            },
        ],
    },

    ##################################################
    # formats: MAF
    ##################################################

    {
        'group'   => "MAF",
        'header'  => "UCSC MAF:",
        'options' => [
            {
                'option'  => "block",
                'type'    => "align::block",
            },
        ],
    },

    ##################################################
    # formats: MULTAL
    ##################################################

    {
        'group'   => "MULTAL",
        'header'  => "MULTAL/MULTAS:",
        'options' => [
            {
                'option'  => "block",
                'type'    => "align::block",
            },
        ],
    },

    ##################################################
    # load user colormaps / groupmaps
    ##################################################

    {
        'group'   => "MAPFILES",
        'header'  => "User defined colormap and consensus group definition:",
        'options' => [
            {
                'option'  => "colorfile",
                'usage'   => "Load more colormaps from file",
                'type'    => "infile",
                'order'   => -1,  #done before -listcolors
                'action'  => sub {
                    my ($self, $pn, $pv, $e) = @_;
                    local *TMP;
                    return  unless $pv;
                    return  unless open(TMP, "< $pv");
                    Bio::MView::Manager::load_colormaps(\*TMP);
                    close TMP;
                },
            },

            {
                'option'  => "groupfile",
                'usage'   => "Load more groupmaps from file",
                'type'    => "infile",
                'order'   => -1,  #done before -listgroups
                'action'  => sub {
                    my ($self, $pn, $pv, $e) = @_;
                    local *TMP;
                    return  unless $pv;
                    return  unless open(TMP, "< $pv");
                    Bio::MView::Manager::load_groupmaps(\*TMP);
                    close TMP;
                },
            },
        ],
    },

    ##################################################
    # info and help
    ##################################################

    {
        'group'   => "INFO",
        'header'  => "More information and help:",
        'options' => [
            {
                'option'  => "help",
                'usage'   => "This help",
                'type'    => "flag::silent",
                'order'   => -1,
            },

            {
                'option'  => "listcolors",
                'usage'   => "Print listing of known colormaps",
                'type'    => "flag::silent",
                'order'   => 0,
            },

            {
                'option'  => "listgroups",
                'usage'   => "Print listing of known consensus groups",
                'type'    => "flag::silent",
                'order'   => 0,
            },

            {
                'option'  => "listcss",
                'usage'   => "Print style sheet",
                'type'    => "flag::silent",
                'order'   => 0,
            },
        ],
    },

];

###########################################################################
1;
