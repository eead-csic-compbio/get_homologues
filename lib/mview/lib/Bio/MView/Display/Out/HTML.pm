# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::MView::Display::Out::HTML;

use vars qw(@ISA);

@ISA = qw(Bio::MView::Display::Out::Text);

######################################################################
# public methods
######################################################################
#overrides
sub process_char {
    my ($self, $c, $caller) = @_;

    my $i = $caller->{'cursor'};

    my $class = $caller->{'r_class'}->[$i];
    if (defined $class) {
        return ("<SPAN CLASS=$class>", $c, "</SPAN>");
    }

    my $color = $caller->{'r_color'}->[$i];
    if (defined $color) {
        return "<SPAN style=\"color:$color\">", $c, "</SPAN>";
    }

    return ($c);
}

#overrides
sub process_segment {
    my ($self, $string) = @_;
    return strip_html_repeats($string);
}

###########################################################################
#overrides
sub render_identifier {
    my ($self, $w, $s, $url) = @_;
    my $rpad = ' ' x ($w - length $s);
    $s = "<A HREF=\"$url\">$s</A>"  if defined $url and $url;
    $self->write($s, $rpad);
}

#overrides
sub render_position {
    my ($self, $just, $w, $s, $isruler, $bold) = @_;
    $s = $self->format_label($just, $w, ($isruler ? $s : ''));
    $s = "<STRONG>$s</STRONG>"  if $bold;
    $self->write($s);
}

#overrides
sub render_sequence {
    my ($self, $s, $bold) = @_;
    $s = "<STRONG>$s</STRONG>"  if $bold;
    $self->write($s);
}

#overrides
sub render_table_begin {
    my $self = shift;
    my $par = shift;

    my $alncolor   = $par->{'alncolor'};
    my $labcolor   = $par->{'labcolor'};
    my $linkcolor  = $par->{'linkcolor'};
    my $alinkcolor = $par->{'alinkcolor'};
    my $vlinkcolor = $par->{'vlinkcolor'};

    #table attrs
    my $s = "style=\"border:0px;";
    if (! $par->{'css1'}) {
        #supported in HTML 4.01:
        $s .= " background-color:$alncolor;"  if defined $alncolor;
        $s .= " color:$labcolor;"             if defined $labcolor;
        $s .= " a:link:$linkcolor;"           if defined $linkcolor;
        $s .= " a:active:$alinkcolor;"        if defined $alinkcolor;
        $s .= " a:visited:$vlinkcolor;"       if defined $vlinkcolor;
    }
    $s .= "\"";
    $self->write("<TABLE $s>\n");
}

#overrides: any undef will be ignored
sub render_text {
    my $self = shift;
    foreach my $s (@_) {
        $self->write($s), next  if defined $s;
        #ignore undef
    }
}

#overrides
sub render_table_end {
    my $self = shift;
    $self->write("</TABLE>\n");
}

#overrides
sub render_tr_pre_begin {
    my $self = shift;
    $self->write("<TR><TD><PRE>\n");
}

#overrides
sub render_tr_pre_end {
    my $self = shift;
    $self->write("</PRE></TD></TR>\n");
}

######################################################################
# private class methods
######################################################################
sub strip_html_repeats {
    my ($list) = @_;

    return $list  if @$list < 4;

    #warn " IN=[", @$list, "]\n";

    #use 1 element lookahead
    my ($limit, $new) = (scalar @$list, []);

    for (my $i=0; $i<$limit; $i++) {

        if ($list->[$i] =~ /^</) {
            #looking at: <tag> or </tag>

            #empty stack
            if (@$new < 1) {
                #warn "NEW(@{[scalar @$new]}) lex[$i] <-- $list->[$i]\n";
                push @$new, $list->[$i];
                #warn "[", @$new, "]\n";
                next;
            }

            #non-empty stack: different tag
            if ($list->[$i] ne $new->[$#$new]) {
                #warn "NEW(@{[scalar @$new]}) lex[$i] <-- $list->[$i]\n";
                push @$new, $list->[$i];
                #warn "[", @$new, "]\n";
                next;
            }

            #non-empty stack: same tag
            #warn "NOP(@{[scalar @$new]}) lex[$i]     $new->[$#$new]\n";
            #warn "[", @$new, "]\n";

        } else {
            #non-tag
            #warn "ARG(@{[scalar @$new]})    [$i] <-- $list->[$i]\n";
            push @$new, $list->[$i];
            #warn "[", @$new, "]\n";
        }
    }

    #warn "SR1=[", @$new, "]\n";

    $new = strip_html_tandem_repeats($new);

    #warn "SR2=[", @$new, "]\n\n";

    return $new;
}

sub strip_html_tandem_repeats {
    my ($list) = @_;

    return $list  if @$list < 6;

    #warn " IN=[", @$list, "]\n";

    #use 1 element lookahead
    my ($limit, $new, @mem) = (scalar @$list, []);

    for (my $i=0; $i<$limit; $i++) {

        #looking at: tag
        if ($list->[$i] =~ /^<[^\/]/) {

            #empty stack
            if (@mem < 1) {
                #warn "NEW(@{[scalar @mem]}) lex[$i] <-- $list->[$i]\n";
                push @mem, $list->[$i];
                #scan until next closing tag
                $i = close_tag($list, $limit, ++$i, \@mem);
                next;
            }

            #non-empty stack: different tag
            if ($list->[$i] ne $mem[0]) {
                #warn "NEW(@{[scalar @mem]}) lex[$i] <-- $list->[$i]\n";
                push @$new, @mem;
                @mem = ();
                redo;
            }

            #non-empty stack: same tag
            #warn "POP(@{[scalar @mem]}) lex[$i] --> $mem[$#mem]\n";
            pop @mem;
            #scan until next closing tag
            $i = close_tag($list, $limit, ++$i, \@mem);
            next;

        } else {
            #non-tag
            #warn "ARG(@{[scalar @mem]}) lex[$i] <-- $list->[$i]\n";
            push @$new, @mem, $list->[$i];
            @mem = ();
        }
    }

    push @$new, @mem    if @mem;

    #warn "SR2=[", @$new, "]\n\n";

    return $new;
}

sub close_tag {
    my ($list, $limit, $i, $mem) = @_;
    my $tag = 0;
    while ($i < $limit) {

        if ($list->[$i] =~ /^<[^\/]/) {
            #inner tag
            $tag++;
            #warn "I  tag[$i] <-- $list->[$i]\n";
            push @$mem, $list->[$i++];
            next;
        }

        if ($list->[$i] =~ /^<\//) {
            if ($tag) {
                #inner /tag
                $tag--;
                #warn "I /tag[$i] <-- $list->[$i]\n";
                push @$mem, $list->[$i++];
                next;
            }
            #outer /tag
            #warn "O /tag[$i] <-- $list->[$i]\n";
            push @$mem, $list->[$i];
            last;
        }

        #datum
        #warn "  data[$i] <-- $list->[$i]\n";
        push @$mem, $list->[$i++];
    }
    return $i;
}

###########################################################################
1;
