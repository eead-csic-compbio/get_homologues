# Copyright (C) 1997-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::MView::Manager;

use Bio::Util::System qw(vmstat);
use Bio::MView::Option::Parameters;  #for $PAR
use Bio::MView::Option::Arguments;
use Bio::MView::Color::ColorMap;
use Bio::MView::Build::Base;
use Bio::MView::Convert;
use Bio::Parse::Stream;

sub new {
    my $type = shift;
    die "${type}::new() missing argument\n"  if @_ < 1;
    my $self = {};
    bless $self, $type;

    $self->{'acount'}  = 0;      #alignments processed
    $self->{'fcount'}  = 0;      #files processed
    $self->{'file'}    = undef;
    $self->{'format'}  = undef;
    $self->{'stream'}  = undef;
    $self->{'class'}   = undef;
    $self->{'build'}   = undef;
    $self->{'display'} = shift;  #display object

    $self;
}

######################################################################
# public class methods
######################################################################
sub check_input_file {
    my $file = shift;
    return Bio::MView::Option::Arguments::check_informat($file, 'file');
}

sub load_colormaps { Bio::MView::Color::ColorMap::load_colormaps(@_) }
sub dump_colormaps { Bio::MView::Color::ColorMap::dump_colormaps(@_) }
sub dump_css       { Bio::MView::Color::ColorMap::dump_css1(@_) }

sub load_groupmaps { Bio::MView::GroupMap::load_groupmaps(@_) }
sub dump_groupmaps { Bio::MView::GroupMap::dump_groupmaps(@_) }

######################################################################
# public methods
######################################################################
#Called with the desired format to be parsed: either a string 'X' naming a
#Parse::Format::X or a hint which will be recognised by that class.
sub parse {
    my ($self, $file, $format) = @_;

    $self->{'file'}   = $file;
    $self->{'format'} = lc $format;
    $self->{'class'}  = load_format_library($format);
    $self->{'stream'} = get_parser_stream($file, $self->{'class'});

    # warn $self->{'file'}, "\n";
    # warn $self->{'format'}, "\n";
    # warn $self->{'class'}, "\n";
    # warn $self->{'stream'}, "\n";

    return 0  unless defined $self->{'stream'};

    my $dis = $self->{'display'};

    while (my $bld = $self->next_build) {

        last  unless defined $bld;  #all done

        $self->{'fcount'}++;  #file count

        $bld->reset;

        my $ac = 0;  #local aln count

        while (defined (my $aln = $bld->next_align)) {

            if ($aln < 1) {
                #warn "parse: empty alignment\n";
                next;
            }

            $self->{'acount'}++;  #global aln count

            $aln->sort_alignment($PAR->get('sort'));

            if ($PAR->get('outfmt') ne 'mview') {
                $self->print_format_conversion($PAR, $bld, $aln);
            } else {
                $self->add_alignment_display($bld, $aln, $dis, ++$ac);

                #display item now?
                $self->{'display'}->render  unless $PAR->get('register');
            }

            $aln = undef;  #gc
            #vmstat("Manager: Align dropped");
        }

        $bld = undef;  #gc
        #vmstat("Manager: Build dropped");
    }

    return 1;
}

sub get_alignment_count { return $_[0]->{'acount'} }

######################################################################
# private class methods
######################################################################
sub load_format_library {
    my $class = "Bio::MView::Build::Format::$_[0]";
    my $library = $class;
    $library =~ s/::/\//g;
    require "$library.pm";
    return $class;
}

sub get_format_parser {
    my $parser = $_[0] . "::parser";
    no strict 'refs';
    $parser = &$parser();
    use strict 'refs';
    return $parser;
}

sub get_parser_stream {
    my ($file, $class) = @_;
    my $parser = get_format_parser($class);
    return new Bio::Parse::Stream($file, $parser);
}

######################################################################
# private methods
######################################################################
#return next entry worth of parse data as in a Bio::MView::Build::Base
#derived object ready for parsing, or undef if no more data.
sub next_build {
    my $self = shift;

    #free last entry and garbage its Bio::MView::Build derived object
    if (defined $self->{'build'}) {
        $self->{'build'}->get_entry->free;
        $self->{'build'} = undef;
    }

    #read the next chunk of data
    my $entry = $self->{'stream'}->get_entry;
    if (! defined $entry) {
        $self->{'stream'}->close;
        return undef;
    }

    #construct a new Bio::MView::Build derived object
    return $self->{'build'} = $self->{'class'}->new($entry);
}

sub gc_flag {
    return 0  if $PAR->get('consensus');
    return 0  if $PAR->get('conservation');
    return 1;
}

#construct a header string describing this alignment
sub header {
    my $self = shift;
    return "Input: $self->{'file'}\n";
}

sub assemble_headers {
    my ($self, $bld, $aln, $acount) = @_;

    my $nfiles  = scalar @{$PAR->get('argvfiles')};
    my $fcount = $self->{'fcount'};
    my @headers = ();
    my $s;

    #first alignment
    if ($acount < 2) {

        #multiple input files: filename on first alignment
        if ($nfiles > 1) {
            $s = $self->header;
            if ($fcount > 1) {
                push @headers, undef; #mark start of new file
            }
            push @headers, $s;
        }

        #build and alignment headers on first alignment
        $s = $bld->header; push @headers, $s  if $s ne '';
        $s = $aln->header; push @headers, $s  if $s ne '';
    }

    #subheaders
    $s = $bld->subheader; push @headers, "\n", $s  if $s ne '';
    push @headers, undef  if @headers;  #mark end of headers

    #warn "(fcount: $fcount, acount: $acount)\n";
    #my @tmp=map {$_=defined $_?$_:'undef';s/\n/\\n/g;"\"$_\""} @{[@headers]};
    #warn "assemble_headers: [@{[join(',', @tmp)]}]\n";

    return \@headers;
}

sub add_alignment_display {
    my ($self, $bld, $aln, $dis, $acount) = @_;

    my $headers = $self->assemble_headers($bld, $aln, $acount);

    #vmstat("display panel constructor");
    my $panel = $dis->new_panel($headers, $aln->init_display);
    #vmstat("display panel constructor DONE");

    my $refobj = $bld->get_ref_row;
    my $refuid = $bld->get_ref_uid;

    #attach a ruler? (may include header text)
    if ($PAR->get('ruler')) {
        my $tmp = $aln->build_ruler($refobj);
        $tmp->append_display($panel);
        #vmstat("ruler added");
    }

    #attach the alignment
    if ($PAR->get('alignment')) {
        $aln->set_color_scheme($refuid);
        $aln->append_display($panel, $self->gc_flag);
        #vmstat("alignment added");
    }

    #attach conservation line?
    if ($PAR->get('conservation')) {
        my $tmp = $aln->build_conservation_row;
        $tmp->append_display($panel);
        #vmstat("conservation added");
    }

    #attach consensus alignments?
    if ($PAR->get('consensus')) {
        my $tmp = $aln->build_consensus_rows;
        $tmp->set_consensus_color_scheme($aln, $refuid);
        $tmp->append_display($panel);
        #vmstat("consensus added");
    }

    #garbage collect if not already done piecemeal
    if (!$self->gc_flag) {
        $aln->do_gc;
        #vmstat("final garbage collect");
    }
}

sub print_format_conversion {
    my ($self, $par, $bld, $aln, $stm) = (@_, \*STDOUT);
    my $conv = new Bio::MView::Convert($bld, $aln, $par->get('moltype'));
    my $outfmt = $par->get('outfmt');
    my $s;
    while (1) {
        $s = $conv->clustal,  last  if $outfmt eq 'clustal';
        $s = $conv->msf,      last  if $outfmt eq 'msf';
        $s = $conv->pearson,  last  if $outfmt eq 'pearson';
        $s = $conv->pir,      last  if $outfmt eq 'pir';
        $s = $conv->plain,    last  if $outfmt eq 'plain';
        $s = $conv->rdb,      last  if $outfmt eq 'rdb';
        last;
    }
    print $stm $$s  if defined $s;
}

###########################################################################
1;
