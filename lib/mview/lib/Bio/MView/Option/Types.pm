# Copyright (C) 2018-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::MView::Option::Types;

use Bio::MView::Option::Arguments;
use Bio::MView::Align::Consensus;
use Bio::MView::Align::Sequence;
use Bio::MView::Color::ColorMap;
use Bio::MView::GroupMap;
use Bio::Util::Regexp;
use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);

push @EXPORT, qw(
    get_default_sequence_colormap
    get_default_consensus_colormap
    get_default_groupmap
);

use vars qw($Types);

######################################################################
sub get_default_sequence_colormap {
    Bio::MView::Color::ColorMap::get_default_sequence_colormap(@_)
}
sub get_default_consensus_colormap {
    Bio::MView::Color::ColorMap::get_default_consensus_colormap(@_)
}
sub get_default_groupmap {
    Bio::MView::GroupMap::get_default_groupmap(@_)
}

######################################################################
$Types = [

    #### REUSABLE TYPES ####

    {
        'type'    => "flag",  #no argument
        'label'   => "",
        'default' => "unset",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return 0  if $ov eq 'unset';  #default
            return 1  if $ov eq 'set';
            return 0  if $ov eq '0';
            return 1  if $ov eq '1';
        },
    },

    {
        'type'    => "flag::silent",  #no argument, empty default
        'label'   => "",
        'default' => "",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return 0  if $ov eq '';  #default
            return 0  if $ov eq '0';
            return 1  if $ov eq '1';
        },
    },

    {
        'type'    => "flag::invert",  #no argument, invert input
        'label'   => "",
        'default' => "set",  #preset
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return 0  if $ov eq 'unset';  #confirm preset
            return 1  if $ov eq 'set';    #confirm preset
            return 1  if $ov eq '0';      #invert input
            return 0  if $ov eq '1';      #invert input
            return 1;  #same as default
        },
    },

    {
        'type'    => "binary",
        'label'   => "on|off",
        'default' => "off",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return 0  if $ov eq 'off' or $ov eq '0';
            return 1  if $ov eq 'on'  or $ov eq '1';
            push(@$e, "bad $on value '$ov', want on|off");
            return $ov;
        },
    },

    {
        'type'    => "char",
        'label'   => "char",
        'default' => "",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            push(@$e, "bad $on value '$ov', want character"),
                unless length($ov) == 1;
            return $ov;
        },
    },

    {
        'type'    => "string",
        'label'   => "string",
        'default' => "",
    },

    {
        'type'    => "color",  #really just a string
        'label'   => "color",
        'default' => "",
    },

    {
        'type'    => "u_int",
        'label'   => "integer",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            push(@$e, "bad $on value '$ov', want integer >= 0"),
                unless $ov =~ /^$RX_Uint$/;
            return $ov;
        },
    },

    {
        'type'    => "u_numeric",
        'label'   => "number",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            push(@$e, "bad $on value '$ov', want numeric >= 0"),
                unless $ov =~ /^$RX_Ureal$/;
            return $ov;
        },
    },

    {
        'type'    => "numeric_gt_0",
        'label'   => "N",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = $self->test_type('u_numeric', $on, $ov, $e);
            my $msg = "bad value '$ov', want numeric > 0";
            if (! defined $pv) { pop @$e; push(@$e, $msg); return undef }
            if ($pv <= 0)      {          push(@$e, $msg); return undef }
            return $pv;
        },
    },

    {
        'type'    => "percentage",
        'label'   => "N",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = $self->test_type('u_numeric', $on, $ov, $e);
            return $pv  if 0 <= $pv and $pv <= 100;
            push(@$e, "bad $on value '$ov', want percentage in 0..100");
            return $ov;
        },
    },

    {
        'type'    => "percentage_ge_50",
        'label'   => "N",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = $self->test_type('percentage', $on, $ov, $e);
            return $pv  if 50 <= $pv and $pv <= 100;
            push(@$e, "bad $on value '$ov', want percentage in 50..100");
            return $ov;
        },
    },

    {
        'type'    => "infile",
        'label'   => "file",
        'default' => "",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            if (defined $ov and $ov ne '') {
                local *TMP;
                open(TMP, "< $ov") or push(@$e, "can't open $on file '$ov'");
                close TMP;
            }
            return $ov;
        },
    },

    {
        'type'    => "string_list",
        'label'   => "str[,str]",
        'default' => "",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return []  unless defined $ov;
            my @tmp = ();
            local $_;
            #warn "type: string_list($on, $ov) in\n";
            foreach (split /[,\s]+/, $ov) {
                next  unless length($_);
                #integer range M..N or M:N
                if (/^($RX_Uint)(?:\.\.|:)($RX_Uint)$/) {
                    if ($2 < $1) {
                        push @tmp, $2..$1;
                    } else {
                        push @tmp, $1..$2;
                    }
                    next;
                }
                #non-range: take whole string
                push @tmp, $_;
            }
            #warn "type: string_list(@tmp) out\n";
            return [ @tmp ];
        },
    },

    {
        'type'    => "u_int_list",
        'label'   => "int[,int]",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return []  unless defined $ov;
            my @tmp = ();
            local $_;
            #warn "type: u_int_list($on, $ov)\n";
            foreach (split /[,\s]+/, $ov) {
                next  unless length($_);
                #positive integer range M..N or M:N
                if (/^($RX_Uint)(?:\.\.|:)($RX_Uint)$/) {
                    if ($2 < $1) {
                        push @tmp, $2..$1;
                    } else {
                        push @tmp, $1..$2;
                    }
                    next;
                }
                #non-range
                if (/^($RX_Uint)$/ and ! /\.\./ and ! /:/) {
                    push @tmp, $1;
                    next;
                }
                push @$e, "bad $on value '$ov', want integer list";
                return [];
            }
            #warn "test: u_int_list(@tmp)\n";
            return [ @tmp ];
        },
    },

    {
        'type'    => "u_numeric_list",
        'label'   => "num[,num]",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return []  unless defined $ov;
            my @tmp = ();
            local $_;
            #warn "type: u_numeric_list($on, $ov)\n";
            foreach (split /[,\s]+/, $ov) {
                next  unless length($_);
                #non-range
                if (/^($RX_Ureal)$/ and ! /\.\./ and ! /:/) {
                    push @tmp, $1;
                    next;
                }
                push @$e, "bad $on value '$ov', want numeric list";
                return [];
            }
            #warn "type: u_numeric_list(@tmp)\n";
            return [ @tmp ];
        },
    },

    {
        'type'    => "integer_range",
        'label'   => "M:N,all",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return []  if $ov eq 'all';
            my @tmp = split(/:/, $ov);
            my $msg = "bad $on range '$ov', want M:N,all";
            if (@tmp != 2) {
                push @$e, $msg; return [];
            }
            $self->test_type('u_int', $on, $tmp[0], $e);
            if (@$e) {
                push @$e, $msg; return [];
            }
            $self->test_type('u_int', $on, $tmp[1], $e);
            if (@$e) {
                push @$e, $msg; return [];
            }
            if ($tmp[0] > $tmp[1]) {  #ensure ascending range
                my $tmp = $tmp[0]; $tmp[0] = $tmp[1]; $tmp[1] = $tmp;
            }
            return [ @tmp ];
        },
    },

    #### MVIEW TYPES ####

    {
        'type'    => "format::in",
        'label'   => "format",
        'values'  => list_informats,
        'default' => "",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_informat($ov);
            if (! defined $pv) {
                push @$e, "input format '$ov' unrecognised";
                push @$e, "valid values are: @{[list_informats]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "format::out",
        'label'   => "format",
        'values'  => list_outformats,
        'default' => "mview",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_outformat($ov);
            if (! defined $pv) {
                push @$e, "output format '$ov' unrecognised";
                push @$e, "valid values are: @{[list_outformats]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "content::width",
        'label'   => "columns",
        'values'  => "N,full",
        'default' => "full",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return 0  if $ov eq "full";
            return 0  if $ov eq "flat";  #backwards compatibility
            my $pv = $self->test_type('u_int', $on, $ov, $e);
            if (@$e) {
                pop @$e;
                push @$e, "width '$ov' unrecognised";
                push @$e, "valid values are: {N,full} where N is an integer";
            }
            return $pv;
        },
    },

    {
        'type'    => "identity::pcid",
        'label'   => "mode",
        'values'  => list_identity_modes,
        'default' => "aligned",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_identity_mode($ov);
            if (! defined $pv) {
                push @$e, "percent identity mode '$ov' unrecognised";
                push @$e, "valid values are: @{[list_identity_modes]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "sort::mode",
        'label'   => "mode",
        'values'  => list_sort_modes,
        'default' => "none",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_sort_mode($ov);
            if (! defined $pv) {
                push @$e, "sort mode '$ov' unrecognised";
                push @$e, "valid values are: @{[list_sort_modes]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "filter::top",
        'label'   => "count",
        'values'  => "N,all",
        'default' => "all",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return 0  if $ov eq "all";
            my $pv = $self->test_type('u_int', $on, $ov, $e);
            if (@$e) {
                pop @$e;
                push @$e, "topn value '$ov' unrecognised";
                push @$e, "valid values are: {N,all} where N is an integer";
            }
            return $pv;
        },
    },

    {
        'type'    => "moltype::type",
        'label'   => "mol",
        'values'  => list_molecule_types,
        'default' => "aa",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_molecule_type($ov);
            if (! defined $pv) {
                push @$e, "molecule type '$ov' unrecognised";
                push @$e, "valid values are: @{[list_molecule_types]}";
                return $ov;
            }
            return $pv;
        },
    },

    {
        'type'    => "alignment::coloring",
        'label'   => "mode",
        'values'  => list_alignment_color_schemes,
        'default' => "none",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_alignment_color_scheme($ov);
            if (! defined $pv) {
                push @$e, "alignment coloring scheme '$ov' unrecognised";
                push @$e, "valid values are: @{[list_alignment_color_schemes]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "alignment::colormap",
        'label'   => "name",
        #'values'  => list_colormap_names,
        'values'  => "see manual",
        'default' => get_default_sequence_colormap,
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_colormap($ov);
            if (! defined $pv) {
                push @$e, "alignment colormap '$ov' unrecognised";
                push @$e, "valid values are: @{[list_colormap_names]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "alignment::groupmap",
        'label'   => "name",
        #'values'  => list_groupmap_names,
        'values'  => "see manual",
        'default' => get_default_groupmap,
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_groupmap($ov);
            if (! defined $pv) {
                push @$e, "alignment groupmap '$ov' unrecognised";
                push @$e, "valid values are: @{[list_groupmap_names]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "alignment::ignore",
        'values'  => list_ignore_classes,
        'label'   => "mode",
        'default' => "none",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_ignore_class($ov);
            if (! defined $pv) {
                push @$e, "ignore class '$ov' unrecognised";
                push @$e, "valid values are: @{[list_ignore_classes]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "consensus::coloring",
        'values'  => list_consensus_color_schemes,
        'label'   => "mode",
        'default' => "none",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_consensus_color_scheme($ov);
            if (! defined $pv) {
                push @$e, "consensus coloring scheme '$ov' unrecognised";
                push @$e, "valid values are: @{[list_consensus_color_schemes]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "consensus::colormap",
        #'values'  => list_colormap_names,
        'values'  => "see manual",
        'label'   => "name",
        'default' => get_default_consensus_colormap,
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_colormap($ov);
            if (! defined $pv) {
                push @$e, "consensus colormap '$ov' unrecognised";
                push @$e, "valid values are: @{[list_colormap_names]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "consensus::groupmap",
        #'values'  => list_groupmap_names,
        'values'  => "see manual",
        'label'   => "name",
        'default' => get_default_groupmap,
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_groupmap($ov);
            if (! defined $pv) {
                push @$e, "consensus groupmap '$ov' unrecognised";
                push @$e, "valid values are: @{[list_groupmap_names]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "consensus::threshold",
        'label'   => "N[,N]",
        'default' => '100,90,80,70',
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = $self->test_type('u_numeric_list', $on, $ov, $e);
            foreach my $v (@$pv) {
                if ($v < 50 or $v > 100) {
                    push @$e, "bad $on value '$ov', want percentage in range 50..100";
                    return undef;
                }
            }
            return $pv;
        },
    },

    {
        'type'    => "consensus::ignore",
        'values'  => list_ignore_classes,
        'label'   => "mode",
        'default' => "none",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_ignore_class($ov);
            if (! defined $pv) {
                push @$e, "consensus ignore class '$ov' unrecognised";
                push @$e, "valid values are: @{[list_ignore_classes]}";
            }
            return  $pv;
        },
    },

    {
        'type'    => "html::html_mode",
        'values'  => list_html_modes,
        'label'   => "mode",
        'default' => "off",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_html_mode($ov);
            if (! defined $pv) {
                push @$e, "html mode '$ov' unrecognised";
                push @$e, "valid values are: @{[list_html_modes]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "html::css_mode",
        'values'  => list_css_modes,
        'label'   => "mode",
        'default' => "off",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_css_mode($ov);
            if (! defined $pv) {
                push @$e, "css mode '$ov' unrecognised";
                push @$e, "valid values are: @{[list_css_modes]}";
            }
            return $pv;
        },
    },

    #### FORMAT SPECIFIC TYPES ####

    {
        'type'    => "blast::hsp_mode",
        'usage'   => "HSP tiling mode",
        'values'  => list_hsp_tiling,
        'label'   => "mode",
        'default' => "ranked",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_hsp_tiling($ov);
            if (! defined $pv) {
                push @$e, "bad hsp mode '$ov'";
                push @$e, "valid values are: @{[list_hsp_tiling]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "blast::maxpval",
        'usage'   => "Ignore hits with p-value greater than N",
        'default' => "unlimited",
        'label'   => "N,unlimited",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return undef  if $ov eq 'unlimited';
            my $pv = $self->test_type('numeric_gt_0', $on, $ov, $e);
            if (! defined $pv) {
                pop @$e;
                push(@$e, "bad maxpval value '$ov', want numeric > 0");
            }
            return $pv;
        },
    },

    {
        'type'    => "blast::minscore",
        'usage'   => "Ignore hits with score less than N",
        'label'   => "N,unlimited",
        'default' => "unlimited",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return undef  if $ov eq 'unlimited';
            my $pv = $self->test_type('numeric_gt_0', $on, $ov, $e);
            if (! defined $pv) {
                pop @$e;
                push(@$e, "bad minscore value '$ov', want numeric > 0");
            }
            return $pv;
        },
    },

    {
        'type'    => "blast::maxeval",
        'usage'   => "Ignore hits with e-value greater than N",
        'label'   => "N,unlimited",
        'default' => "unlimited",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return undef  if $ov eq 'unlimited';
            my $pv = $self->test_type('numeric_gt_0', $on, $ov, $e);
            if (! defined $pv) {
                pop @$e;
                push(@$e, "bad maxeval value '$ov', want numeric > 0");
            }
            return $pv;
        },
    },

    {
        'type'    => "blast::minbits",
        'usage'   => "Ignore hits with bits less than N",
        'label'   => "N,unlimited",
        'default' => "unlimited",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return undef  if $ov eq 'unlimited';
            my $pv = $self->test_type('numeric_gt_0', $on, $ov, $e);
            if (! defined $pv) {
                pop @$e;
                push(@$e, "bad minbits value '$ov', want numeric > 0");
            }
            return $pv;
        },
    },

    {
        'type'    => "blast::keepinserts",
        'usage'   => "Keep hit sequence insertions in unaligned output",
        'label'   => "on|off",
        'default' => "off",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return $self->test_type('binary', $on, $ov, $e);
        },
    },

    {
        'type'    => "psiblast::cycle",
        'usage'   => "Process the N'th cycle of a multipass search",
        'values'  => list_cycle_types,
        'label'   => "cycles",
        'default' => "last",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_cycle_type($ov);
            if (! defined $pv) {
                push @$e, "bad cycle value in '$ov'";
                push @$e, "valid values are: @{[list_cycle_types]}";
            }
            return $pv;
        },
    },


    {
        'type'    => "search::strand",
        'usage'   => "Report only these query strand orientations",
        'values'  => list_strand_types,
        'label'   => "strands",
        'default' => "both",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_strand_type($ov);
            if (! defined $pv) {
                push @$e, "bad strand orientation in '$ov'";
                push @$e, "valid values are: @{[list_strand_types]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "align::block",
        'usage'   => "Report only these blocks",
        'values'  => list_block_values,
        'label'   => "blocks",
        'default' => "all",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_block_value($ov);
            if (! defined $pv) {
                push @$e, "bad block value in '$ov'";
                push @$e, "valid values are: @{[list_block_values]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "hssp::chain",
        'usage'   => "Report only these chain names/numbers",
        'values'  => list_chain_values,
        'label'   => "chains",
        'default' => "all",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            my $pv = check_chain_value($ov);
            if (! defined $pv) {
                push @$e, "bad chain in '$ov'";
                push @$e, "valid values are: @{[list_chain_values]}";
            }
            return $pv;
        },
    },

    {
        'type'    => "fasta::minopt",
        'usage'   => "Ignore hits with opt score less than N",
        'label'   => "N,unlimited",
        'default' => "unlimited",
        'test'    => sub {
            my ($self, $on, $ov, $e) = @_;
            return undef  if $ov eq 'unlimited';
            my $pv = $self->test_type('numeric_gt_0', $on, $ov, $e);
            if (! defined $pv) {
                pop @$e;
                push(@$e, "bad minscore value '$ov', want numeric > 0");
            }
            return $pv;
        },
    },

];

###########################################################################
1;
