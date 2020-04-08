# Copyright (C) 1997-2018 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

######################################################################
package Bio::MView::Color::ColorMap;

use Bio::MView::Color::Palette;
use Exporter;

use vars qw(@ISA @EXPORT $COLORMAP);

@ISA = qw(Exporter);

@EXPORT = qw($COLORMAP);

use vars qw($Colour_Black $Colour_White $Colour_DarkGray $Colour_LightGray
            $Colour_Cream $Colour_Brown $Colour_Comment);

$Colour_Black      = '#000000';
$Colour_White      = '#FFFFFF';
$Colour_DarkGray   = '#666666';
$Colour_LightGray  = '#999999';
$Colour_Cream      = '#FFFFCC';
$Colour_Brown      = '#aa6666';
$Colour_Comment    = $Colour_Brown;

my $Palette  = new Bio::MView::Color::Palette;  #static maps names to rgb
my $ColorMap = {};   #static hash of colormaps
my $MapText  = '';   #used as special index
my $Wildcard = '.';  #key for default colouring

my %Default_CSS1_Colors = (
    'alncolor' => $Colour_White,
    'labcolor' => $Colour_Black,
    'symcolor' => $Colour_Black,
    'gapcolor' => $Colour_DarkGray,
);

my $Default_PRO_Aln_ColorMap = 'P1';
my $Default_DNA_Aln_ColorMap = 'D1';
my $Default_PRO_Con_ColorMap = 'PC1';
my $Default_DNA_Con_ColorMap = 'DC1';

$COLORMAP = new Bio::MView::Color::ColorMap;  #unique global instance

sub new {
    my $type = shift;
    if (defined $COLORMAP) {
        die "Bio::MView::Color::ColorMap: instance already exists\n";
    }
    my $self = {};
    bless $self, $type;

    $self->{'map'}     = $ColorMap;
    $self->{'palette'} = $Palette;

    return $COLORMAP = $self;
}

######################################################################
# public class methods
######################################################################
sub colormap_names      { return keys %$ColorMap }

sub list_colormap_names { return join(",", sort keys %$ColorMap) }

sub check_colormap {
    my $map = uc shift;
    return $map  if exists $ColorMap->{$map};   #named colormap
    $map = lc $map;
    return $map  if $Palette->has_color($map);  #predefined colour
    return undef;
}

sub get_default_sequence_colormap {
    return $Default_PRO_Aln_ColorMap  if !defined $_[0];
    return $Default_PRO_Aln_ColorMap  if $_[0] eq 'aa';
    return $Default_DNA_Aln_ColorMap;
}

sub get_default_consensus_colormap {
    return $Default_PRO_Con_ColorMap  if !defined $_[0];
    return $Default_PRO_Con_ColorMap  if $_[0] eq 'aa';
    return $Default_DNA_Con_ColorMap;
}

#load colormap data from stream
sub load_colormaps {
    my ($stream, $override) = (@_, 1);  #allow replacement maps
    my ($state, $map, $de) = ('start');
    local $_;
    while (<$stream>) {
        my $mapignore;

        #comments, blank lines
        if (/^\s*\#/ or /^\s*$/) {
            next  if $state ne 'map';
            $de .= $_;
            next;
        }

        chomp;

        #map [name]
        if (/^\s*\[\s*(\S+)\s*\]/) {
            $state = 'map';
            $map = uc $1;
            if (exists $ColorMap->{$map} and !$override) {
                $mapignore = 1;  #just for duration of this map
            } else {
                $mapignore = 0;
            }
            $de = '';
            next;
        }

        #palette data: colorname: RGBcode
        #colorname MUST begin with a letter
        if (/^\s*colou?r\s+([a-z][-a-z0-9_]*)\s*:\s*(\#[0123456789ABCDEF]{6})/i) {
            $state = 'palette';
            my ($color, $rgb) = (lc $1, lc $2);
            if (! $Palette->has_color($color)) {
                #warn "palette: $color => $rgb\n";
                $Palette->insert($color => $rgb);
            }
            next;
        }

        die "load_colormaps: colormap name undefined\n"  unless defined $map;

        next  if $mapignore;    #forget it if we're not allowing overrides

        #save map description?
        $ColorMap->{$map}->{$MapText} = $de  if $state eq 'map';

        #symbol[symbol] {->|=>} [palette_colorname|#RGB] any other text...
        if (/^\s*(\S)(\S)?\s*(->|=>)\s*(\S+)(?:\s+(.*))?/i) {
            $state = 'symbol';

            my ($c1, $c2, $seethru, $color, $de) =
                ($1, $2, ($3 eq '->' ? 'T' : 'S'),
                 lc $4, (defined $5 ? $5 : ''));

            #only allow new colors in form of RGB codes
            if (! $Palette->has_color($color)) {
                if ($color =~ /\#[0123456789ABCDEF]{6}/i) {
                    #warn "new color: $color (transparency=$seethru)\n";
                    $Palette->insert($color => lc $color);
                } else {
                    die "load_colormaps: undefined color in line '$_'\n";
                }
            }

            #any explicit coloring
            $ColorMap->{$map}->{$c1} =
                [ $Palette->get_index($color), $seethru, $de ];

            $ColorMap->{$map}->{$c2} =
                [ $Palette->get_index($color), $seethru, $de ]
                if defined $c2;

            next;
        }

        #default
        die "load_colormaps: bad format in line '$_'\n";
    }

    close $stream;
}

#return a descriptive listing of all known colormaps
sub dump_colormaps {
    my $html = shift;

    my ($s, $c0, $c1, $c2) = ('', '', '', '');

    ($c0, $c1, $c2) = (
        "<SPAN style=\"color:$Colour_Black\">",
        "<SPAN style=\"color:$Colour_Comment\">",
        "</SPAN>")  if $html;

    sub palette_row {
        my ($f1, $f2, $col, $rgb) = @_;
        return sprintf("color %s%-20s%s : %-7s\n", $f1, $col, $f2, $rgb);
    }

    sub colormap_row {  #binds: $c1, $c2
        my ($p, $key, $f1, $f2, $sym, $col) = @_;
        return sprintf("%s%-7s%s  %s  %s%-20s%s $c1%s$c2\n",
                       $f1, $sym, $f2, ($p->{$key}->[1] eq 'T' ? '->' : '=>'),
                       $f1, $col, $f2, $p->{$key}->[2]);
    }

    $s .= "$c1#Palette:\n\n#color                     : #RGB$c2\n";

    foreach my $color ($Palette->color_names) {
        my ($rgb, $f1, $f2) = ($Palette->get_rgb($color), '', '');
        ($f1, $f2) = ("<SPAN style=\"color:$rgb\">", "</SPAN>")  if $html;
        $s .= palette_row($f1, $f2, $color, $rgb);
    }

    $s .= "\n\n$c1#Colormap listing - suitable for reloading.\n";
    $s .= "#Character matching is case-sensitive.$c2\n\n";

    @_ = keys %$ColorMap  unless @_;

    foreach my $map (sort @_) {
        my %p = %{$ColorMap->{$map}};  #copy colormap structure

        $s .= "$c0\[$map]$c2\n";
        $s .= "$c1$ColorMap->{$map}->{$MapText}";
        $s .= "#symbols =>  color                #comment$c2\n";

        foreach my $sym (sort keys %p) {

            next  if $sym eq $MapText;
            next  unless exists $p{$sym} and defined $p{$sym};

            my $num = $p{$sym}->[0];
            my $col = $Palette->get_color($num);
            my $rgb = $Palette->get_rgb($col);

            #warn "CN: $num,$col,$rgb\n";

            my ($f1, $f2) = ('', '');
            ($f1, $f2) = ("<SPAN style=\"color:$rgb\">", "</SPAN>")  if $html;

            my ($usym, $lsym) = (uc $sym, lc $sym);

            #lower and upper case equal: not alphabetic
            if ($lsym eq $usym) {
                $s .= colormap_row(\%p, $sym, $f1, $f2, $sym, $col);
                next;
            }

            #lower and upper case distinct: alphabetic, both defined
            if (exists $p{$lsym} and exists $p{$usym}) {

                if ($p{$lsym}->[0] eq $p{$usym}->[0] and
                    $p{$lsym}->[1] eq $p{$usym}->[1]) {

                    #common definition: merge them
                    $s .= colormap_row(\%p, $sym, $f1, $f2, $usym.$lsym, $col);

                    $p{$usym} = $p{$lsym} = undef;  #seen them: forget

                } else {  #different definitions; write this one
                    $s .= colormap_row(\%p, $sym, $f1, $f2, $sym, $col);
                }
                next;
            }

            #default: alphabetic, single symbol
            $s .= colormap_row(\%p, $sym, $f1, $f2, $sym, $col);
        }
        $s .= "\n";
    }

    return $s;
}

sub dump_css1 {
    my ($par) = (@_, {});

    my $colors = { %Default_CSS1_Colors };

    $colors->{'alncolor'} = $par->{'alncolor'}  if exists $par->{'alncolor'};
    $colors->{'labcolor'} = $par->{'labcolor'}  if exists $par->{'labcolor'};
    $colors->{'symcolor'} = $par->{'symcolor'}  if exists $par->{'symcolor'};
    $colors->{'gapcolor'} = $par->{'gapcolor'}  if exists $par->{'gapcolor'};

    #warn "bg=$colors->{'alncolor'} fg=$colors->{'symcolor'}\n";

    my $s = "TD {font-family:Fixed,Courier,monospace; background-color:$colors->{'alncolor'}; color:$colors->{'labcolor'}}\n";

    for (my $i=0; $i < $Palette->size; $i++) {

        my $color = $Palette->get_color($i);
        my $rgb   = $Palette->get_rgb_at_index($i);

        #SOLID: coloured background/monochrome foreground: flip foreground
        #between black/white depending on brightest green RGB component of
        #background.
        my $bg = hex(substr($rgb, 1));
        my $fg = (($bg >> 8) & 255) > 200 ? $Colour_Black : $Colour_White;

        #solid: background + foreground
        $s .= "SPAN.S$i \{background-color:$rgb;color:$fg} ";

        #transparent: no background + coloured foreground
        $s .= "SPAN.T$i \{color:$rgb} ";

        #comment
        $s .= "/* $color */\n";
    }

    return $s;
}

######################################################################
# public methods
######################################################################
sub has_colormap {
    my ($self, $map) = @_;
    return 1  if exists $self->{'map'}->{$map};
    return 0;
}

sub has_symbol_color {
    my ($self, $map, $c) = @_;
    return 0  unless exists $self->{'map'}->{$map};
    return 0  unless exists $self->{'map'}->{$map}->{$c};
    return 1;
}

sub has_wildcard_color {
    my ($self, $map) = @_;
    return 0  unless exists $self->{'map'}->{$map}->{$Wildcard};
    return 0  unless exists $self->{'map'}->{$map}->{$Wildcard};
    return 1;
}

sub has_palette_color {
    my ($self, $map) = @_;
    return 1  if $self->{'palette'}->has_color($map);
    return 0;
}

sub get_symbol_color {
    my ($self, $map, $c) = @_;
    my $trans = $self->{'map'}->{$map}->{$c}->[1];
    my $index = $self->{'map'}->{$map}->{$c}->[0];
    my $rgb = $self->{'palette'}->get_rgb_at_index($index);
    return ($rgb, $trans, $index);
}

sub get_wildcard_color {
    my ($self, $map) = @_;
    my $trans = $self->{'map'}->{$map}->{$Wildcard}->[1];
    my $index = $self->{'map'}->{$map}->{$Wildcard}->[0];
    my $rgb = $self->{'palette'}->get_rgb_at_index($index);
    return ($rgb, $trans, $index);
}

sub get_palette_color {
    my ($self, $map, $trans) = (@_, 'S');  #default
    my $index = $self->{'palette'}->get_index($map);
    my $rgb = $self->{'palette'}->get_rgb_at_index($index);
    return ($rgb, $trans, $index);
}

sub get_colormap_length {
    my ($self, $map) = @_;
    return 1  unless exists $self->{'map'}->{$map};
    my $len = scalar keys %{$self->{'map'}->{$map}};
    $len--  if defined $self->{'map'}->{$map}->{''};
    if ($len % 2 != 0) {
        die "get_colormap_length: divide by two error\n"
    }
    return $len / 2;
}

######################################################################
eval { load_colormaps(\*DATA) }; warn($@), exit 1  if $@;

1;

__DATA__

######################################################################
# Original colours with MView
# Netscape 216 (supposedly) cross-platform colours
# (other colours may not resolve well on all platforms)
######################################################################

#Primaries/secondaries (cross-platform)
color  black               :  #000000
color  white               :  #ffffff
color  red                 :  #ff0000
color  green               :  #00ff00
color  blue                :  #0000ff
color  cyan                :  #00ffff
color  magenta             :  #ff00ff
color  yellow              :  #ffff00

#Miscellaneous (cross-platform)
color  purple              :  #6600cc
color  dull-blue           :  #0099ff
color  dark-green-blue     :  #33cccc
color  medium-green-blue   :  #00ffcc
color  bright-blue         :  #0033ff
color  dark-green          :  #009900
color  bright-green        :  #33cc00
color  orange              :  #ff3333
color  orange-brown        :  #cc6600
color  brown               :  #8B4513
color  bright-red          :  #cc0000
color  light-gray          :  #999999
color  dark-gray           :  #666666

#Grey levels (some cross-platform)
color  gray0               :  #ffffff
color  gray1               :  #eeeeee
color  gray2               :  #dddddd
color  gray3               :  #cccccc
color  gray4               :  #bbbbbb
color  gray5               :  #aaaaaa
color  gray6               :  #999999
color  gray7               :  #888888
color  gray8               :  #777777
color  gray9               :  #666666
color  gray10              :  #555555
color  gray11              :  #444444
color  gray12              :  #333333
color  gray13              :  #222222
color  gray14              :  #111111
color  gray15              :  #000000

###########################################################################
# CLUSTALX screen colours (protein + nucleotide)
###########################################################################
color  clustal-red         :  #e53319
color  clustal-blue        :  #197fe5
color  clustal-green       :  #19cc19
color  clustal-cyan        :  #19b2b2
color  clustal-pink        :  #e57f7f
color  clustal-magenta     :  #cc4ccc
color  clustal-yellow      :  #cccc00
color  clustal-orange      :  #e5994c
color  clustal-white       :  #ffffff
color  clustal-light-gray  :  #999999
color  clustal-dark-gray   :  #666666
color  clustal-black       :  #000000

###########################################################################
# CLUSTALX printing colours (protein + nucleotide)
###########################################################################
color  clustal-white-print      :  #ffffff
color  clustal-yellow-print     :  #ffff00
color  clustal-violet-print     :  #6619e5
color  clustal-red-print        :  #e57f66
color  clustal-blue-print       :  #66e5e5
color  clustal-purple-print     :  #b299e5
color  clustal-black-print      :  #000000
color  clustal-pink-print       :  #e57f7f #overrides #cc4ccc in colprint.xml
color  clustal-cyan-print       :  #19b2b2
color  clustal-magenta-print    :  #cc4ccc
color  clustal-orange-print     :  #e5994c #overrides #e5b24c in colprint.xml
color  clustal-green-print      :  #19cc19
color  clustal-light-gray-print :  #999999
color  clustal-dark-gray-print  :  #99b2b2

###########################################################################
# Frederico Nardi's suggested colours for CLUSTAL
###########################################################################
color  nardi-red           :  #ff1111
color  nardi-blue          :  #1155ff
color  nardi-green         :  #11dd11
color  nardi-cyan          :  #11ffff
color  nardi-yellow        :  #ffff11
color  nardi-orange        :  #ff7f11
color  nardi-pink          :  #ff11ff
color  nardi-purple        :  #6611cc
color  nardi-dull-blue     :  #197fe5
color  nardi-light-gray    :  #999999
color  nardi-dark-gray     :  #666666

###########################################################################
# Kuang Lin's colour neural net derived scheme
###########################################################################
color  lin-A               :  #90fe23
color  lin-R               :  #fe5e2d
color  lin-N               :  #2e3d2d
color  lin-D               :  #00903b
color  lin-C               :  #004baa
color  lin-Q               :  #864b00
color  lin-E               :  #3fa201
color  lin-G               :  #10fe68
color  lin-H               :  #b2063b
color  lin-I               :  #04ced9
color  lin-L               :  #4972fe
color  lin-K               :  #c4a100
color  lin-M               :  #2a84dd
color  lin-F               :  #a60ade
color  lin-P               :  #fe61fe
color  lin-S               :  #f7e847
color  lin-T               :  #fefeb3
color  lin-W               :  #4a007f
color  lin-Y               :  #e903a8
color  lin-V               :  #5bfdfd

###########################################################################
# Limited set of MView search block colours
###########################################################################
color  find-A              :  #90fe23
color  find-B              :  #fe5e2d
color  find-C              :  #2e3d2d
color  find-D              :  #00903b
color  find-E              :  #004baa
color  find-F              :  #864b00
color  find-G              :  #3fa201
color  find-H              :  #10fe68
color  find-I              :  #b2063b
color  find-J              :  #04ced9
color  find-K              :  #4972fe
color  find-L              :  #c4a100
color  find-M              :  #2a84dd
color  find-N              :  #a60ade
color  find-O              :  #fe61fe
color  find-P              :  #f7e847
color  find-Q              :  #fefeb3
color  find-R              :  #4a007f
color  find-S              :  #e903a8
color  find-T              :  #5bfdfd

###########################################################################
# Web 4.01 specification:
# https://en.wikipedia.org/wiki/Web_colors
###########################################################################
color  White    : #FFFFFF
color  Silver   : #C0C0C0
color  Gray     : #808080
color  Black    : #000000
color  Red      : #FF0000
color  Maroon   : #800000
color  Yellow   : #FFFF00
color  Olive    : #808000
color  Lime     : #00FF00
color  Green    : #008000
color  Aqua     : #00FFFF
color  Teal     : #008080
color  Blue     : #0000FF
color  Navy     : #000080
color  Fuchsia  : #FF00FF
color  Purple   : #800080

###########################################################################
# X11 colours
# https://en.wikipedia.org/wiki/Web_colors
###########################################################################

#Pink colors
color  Pink                 : #FFC0CB
color  LightPink            : #FFB6C1
color  HotPink              : #FF69B4
color  DeepPink             : #FF1493
color  PaleVioletRed        : #DB7093
color  MediumVioletRed      : #C71585

#Red colors
color  LightSalmon          : #FFA07A
color  Salmon               : #FA8072
color  DarkSalmon           : #E9967A
color  LightCoral           : #F08080
color  IndianRed            : #CD5C5C
color  Crimson              : #DC143C
color  FireBrick            : #B22222
color  DarkRed              : #8B0000
color  Red                  : #FF0000

#Orange colors
color  OrangeRed            : #FF4500
color  Tomato               : #FF6347
color  Coral                : #FF7F50
color  DarkOrange           : #FF8C00
color  Orange               : #FFA500

#Yellow colors
color  Yellow               : #FFFF00
color  LightYellow          : #FFFFE0
color  LemonChiffon         : #FFFACD
color  LightGoldenrodYellow : #FAFAD2
color  PapayaWhip           : #FFEFD5
color  Moccasin             : #FFE4B5
color  PeachPuff            : #FFDAB9
color  PaleGoldenrod        : #EEE8AA
color  Khaki                : #F0E68C
color  DarkKhaki            : #BDB76B
color  Gold                 : #FFD700

#Brown colors
color  Cornsilk             : #FFF8DC
color  BlanchedAlmond       : #FFEBCD
color  Bisque               : #FFE4C4
color  NavajoWhite          : #FFDEAD
color  Wheat                : #F5DEB3
color  BurlyWood            : #DEB887
color  Tan                  : #D2B48C
color  RosyBrown            : #BC8F8F
color  SandyBrown           : #F4A460
color  Goldenrod            : #DAA520
color  DarkGoldenrod        : #B8860B
color  Peru                 : #CD853F
color  Chocolate            : #D2691E
color  SaddleBrown          : #8B4513
color  Sienna               : #A0522D
color  Brown                : #A52A2A
color  Maroon               : #800000

#Green colors
color  DarkOliveGreen       : #556B2F
color  Olive                : #808000
color  OliveDrab            : #6B8E23
color  YellowGreen          : #9ACD32
color  LimeGreen            : #32CD32
color  Lime                 : #00FF00
color  LawnGreen            : #7CFC00
color  Chartreuse           : #7FFF00
color  GreenYellow          : #ADFF2F
color  SpringGreen          : #00FF7F
color  MediumSpringGreen    : #00FA9A
color  LightGreen           : #90EE90
color  PaleGreen            : #98FB98
color  DarkSeaGreen         : #8FBC8F
color  MediumAquamarine     : #66CDAA
color  MediumSeaGreen       : #3CB371
color  SeaGreen             : #2E8B57
color  ForestGreen          : #228B22
color  Green                : #008000
color  DarkGreen            : #006400

#Cyan colors
color  Aqua                 : #00FFFF
color  Cyan                 : #00FFFF
color  LightCyan            : #E0FFFF
color  PaleTurquoise        : #AFEEEE
color  Aquamarine           : #7FFFD4
color  Turquoise            : #40E0D0
color  MediumTurquoise      : #48D1CC
color  DarkTurquoise        : #00CED1
color  LightSeaGreen        : #20B2AA
color  CadetBlue            : #5F9EA0
color  DarkCyan             : #008B8B
color  Teal                 : #008080

#Blue colors
color  LightSteelBlue       : #B0C4DE
color  PowderBlue           : #B0E0E6
color  LightBlue            : #ADD8E6
color  SkyBlue              : #87CEEB
color  LightSkyBlue         : #87CEFA
color  DeepSkyBlue          : #00BFFF
color  DodgerBlue           : #1E90FF
color  CornflowerBlue       : #6495ED
color  SteelBlue            : #4682B4
color  RoyalBlue            : #4169E1
color  Blue                 : #0000FF
color  MediumBlue           : #0000CD
color  DarkBlue             : #00008B
color  Navy                 : #000080
color  MidnightBlue         : #191970

#Purple, violet, and magenta colors
color  Lavender             : #E6E6FA
color  Thistle              : #D8BFD8
color  Plum                 : #DDA0DD
color  Violet               : #EE82EE
color  Orchid               : #DA70D6
color  Fuchsia              : #FF00FF
color  Magenta              : #FF00FF
color  MediumOrchid         : #BA55D3
color  MediumPurple         : #9370DB
color  BlueViolet           : #8A2BE2
color  DarkViolet           : #9400D3
color  DarkOrchid           : #9932CC
color  DarkMagenta          : #8B008B
color  Purple               : #800080
color  Indigo               : #4B0082
color  DarkSlateBlue        : #483D8B
color  SlateBlue            : #6A5ACD
color  MediumSlateBlue      : #7B68EE

#White colors
color  White                : #FFFFFF
color  Snow                 : #FFFAFA
color  Honeydew             : #F0FFF0
color  MintCream            : #F5FFFA
color  Azure                : #F0FFFF
color  AliceBlue            : #F0F8FF
color  GhostWhite           : #F8F8FF
color  WhiteSmoke           : #F5F5F5
color  Seashell             : #FFF5EE
color  Beige                : #F5F5DC
color  OldLace              : #FDF5E6
color  FloralWhite          : #FFFAF0
color  Ivory                : #FFFFF0
color  AntiqueWhite         : #FAEBD7
color  Linen                : #FAF0E6
color  LavenderBlush        : #FFF0F5
color  MistyRose            : #FFE4E1

#Gray and black colors
color  Gainsboro            : #DCDCDC
color  LightGray            : #D3D3D3
color  Silver               : #C0C0C0
color  DarkGray             : #A9A9A9
color  Gray                 : #808080
color  DimGray              : #696969
color  LightSlateGray       : #778899
color  SlateGray            : #708090
color  DarkSlateGray        : #2F4F4F
color  Black                : #000000

######################################################################
# Colourmaps
#
# symbol -> colour (RGB hex or colorname)  [#comment]
######################################################################

[P1]
#Protein: highlight amino acid physicochemical properties
.   ->  dark-gray            #wildcard/mismatch
Gg  =>  bright-green         #hydrophobic
Aa  =>  bright-green         #hydrophobic
Ii  =>  bright-green         #hydrophobic
Vv  =>  bright-green         #hydrophobic
Ll  =>  bright-green         #hydrophobic
Mm  =>  bright-green         #hydrophobic
Ff  =>  dark-green           #large hydrophobic
Yy  =>  dark-green           #large hydrophobic
Ww  =>  dark-green           #large hydrophobic
Hh  =>  dark-green           #large hydrophobic
Cc  =>  yellow               #cysteine
Pp  =>  bright-green         #hydrophobic
Kk  =>  bright-red           #positive charge
Rr  =>  bright-red           #positive charge
Dd  =>  bright-blue          #negative charge
Ee  =>  bright-blue          #negative charge
Qq  =>  purple               #polar
Nn  =>  purple               #polar
Ss  =>  dull-blue            #small alcohol
Tt  =>  dull-blue            #small alcohol
Bb  =>  dark-gray            #D or N
Zz  =>  dark-gray            #E or Q
Xx  ->  dark-gray            #unknown
?   ->  light-gray           #unknown
*   =>  black                #stop

[GPCR]
#Protein: GPCRdb color scheme for Gert Vriend
.   ->  dark-gray            #wildcard/mismatch
Gg  =>  orange-brown         #backbone change
Pp  =>  orange-brown         #backbone change
Aa  =>  bright-green         #hydrophobic
Ii  =>  bright-green         #hydrophobic
Vv  =>  bright-green         #hydrophobic
Ll  =>  bright-green         #hydrophobic
Mm  =>  bright-green         #hydrophobic
Cc  =>  yellow               #cysteine
Qq  =>  bright-red           #positive charge
Ee  =>  bright-red           #positive charge
Nn  =>  bright-red           #positive charge
Dd  =>  bright-red           #positive charge
Bb  =>  bright-red           #D or N
Zz  =>  bright-red           #E or Q
Hh  =>  bright-blue          #negative charge
Kk  =>  bright-blue          #negative charge
Rr  =>  bright-blue          #negative charge
Ss  =>  dark-green-blue      #small alcohol
Tt  =>  dark-green-blue      #small alcohol
Yy  =>  medium-green-blue    #large hydrophobic
Ff  =>  cyan                 #large hydrophobic
Ww  =>  cyan                 #large hydrophobic
Xx  ->  dark-gray            #unknown
?   ->  light-gray           #unknown
*   =>  black                #stop

[CYS]
#Protein: highlight cysteines
.   ->  dark-gray            #wildcard/mismatch
Cc  =>  yellow               #cysteine
Xx  ->  dark-gray            #unknown
?   ->  light-gray           #unknown

[CHARGE]
#Protein: highlight charged amino acids
.   ->  black                #wildcard/mismatch
Kk  =>  bright-red           #positive charge
Rr  =>  bright-red           #positive charge
Dd  =>  bright-blue          #negative charge
Ee  =>  bright-blue          #negative charge
Bb  =>  dark-gray            #D or N
Zz  =>  dark-gray            #E or Q
Xx  =>  dark-gray            #unknown
?   =>  light-gray           #unknown

[POLAR1]
#Protein: highlight charged and polar amino acids
.   ->  dark-gray            #wildcard/mismatch
Kk  =>  bright-red           #positive charge
Rr  =>  bright-red           #positive charge
Dd  =>  bright-blue          #negative charge
Ee  =>  bright-blue          #negative charge
Qq  =>  purple               #charged/polar
Nn  =>  purple               #charged/polar
Ss  =>  purple               #charged/polar
Tt  =>  purple               #charged/polar
Hh  =>  purple               #charged/polar
Bb  =>  purple               #D or N
Zz  =>  purple               #E or Q
Xx  ->  dark-gray            #unknown
?   ->  light-gray           #unknown

[D1]
#DNA: highlight nucleotide types
.   ->  dark-gray            #wildcard/mismatch
Aa  =>  bright-blue          #adenine
Cc  =>  dull-blue            #cytosine
Gg  =>  bright-blue          #guanine
Tt  =>  dull-blue            #thymine
Uu  =>  dull-blue            #uracil
Mm  =>  dark-gray            #amino:      A or C
Rr  =>  dark-gray            #purine:     A or G
Ww  =>  dark-gray            #weak:       A or T
Ss  =>  dark-gray            #strong:     C or G
Yy  =>  dark-gray            #pyrimidine: C or T
Kk  =>  dark-gray            #keto:       G or T
Vv  =>  dark-gray            #not T:      A or C or G
Hh  =>  dark-gray            #not G:      A or C or T
Dd  =>  dark-gray            #not C:      A or G or T
Bb  =>  dark-gray            #not A:      C or G or T
Nn  =>  dark-gray            #any:        A or C or G or T
Xx  ->  dark-gray            #unknown
?   ->  light-gray           #unknown

[D2]
#DNA: highlight match versus mismatch under consensus coloring schemes
.   =>  bright-red           #wildcard/mismatch
Aa  =>  bright-blue          #match
Cc  =>  bright-blue          #match
Gg  =>  bright-blue          #match
Tt  =>  bright-blue          #match
Uu  =>  bright-blue          #match
Mm  =>  bright-blue          #match
Rr  =>  bright-blue          #match
Ww  =>  bright-blue          #match
Ss  =>  bright-blue          #match
Yy  =>  bright-blue          #match
Kk  =>  bright-blue          #match
Vv  =>  bright-blue          #match
Hh  =>  bright-blue          #match
Dd  =>  bright-blue          #match
Bb  =>  bright-blue          #match
Nn  =>  bright-blue          #match
Xx  =>  light-gray           #unknown
?   ->  light-gray           #unknown

[CLUSTAL_NUC]
#CLUSTAL-derived colours for nucleotides
.   =>  clustal-dark-gray    #wildcard/mismatch
Aa  =>  clustal-red          #match
Cc  =>  clustal-blue         #match
Gg  =>  clustal-orange       #match
Tt  =>  clustal-green        #match
Uu  =>  clustal-green        #match
Mm  =>  clustal-dark-gray    #match
Rr  =>  clustal-dark-gray    #match
Ww  =>  clustal-dark-gray    #match
Ss  =>  clustal-dark-gray    #match
Yy  =>  clustal-dark-gray    #match
Kk  =>  clustal-dark-gray    #match
Vv  =>  clustal-dark-gray    #match
Hh  =>  clustal-dark-gray    #match
Dd  =>  clustal-dark-gray    #match
Bb  =>  clustal-dark-gray    #match
Nn  =>  clustal-dark-gray    #match
Xx  =>  clustal-light-gray   #unknown
?   ->  clustal-light-gray   #unknown

[PC1]
#Protein consensus: highlight equivalence class
.  ->  dark-gray             #unconserved
a  ->  dark-green            #aromatic
l  ->  bright-green          #aliphatic
h  ->  bright-green          #hydrophobic
+  ->  bright-red            #positive charge
-  ->  bright-blue           #negative charge
c  ->  purple                #charged
p  ->  dull-blue             #polar
o  ->  dull-blue             #alcohol
u  ->  bright-green          #tiny
s  ->  bright-green          #small
t  ->  bright-green          #turnlike

[DC1]
#DNA consensus: highlight ring type
.  ->  dark-gray             #unconserved
r  ->  purple                #purine
y  ->  orange                #pyrimidine

#+ Frederico Nardi
[NARDI]
#Protein: highlight amino acid physicochemical properties
.        =>  nardi-dark-gray    #wildcard/mismatch
Aa       =>  nardi-dull-blue    #hydrophobic
Bb       =>  nardi-pink         #D or N
Cc       =>  nardi-yellow       #cysteine
Dd       =>  nardi-pink         #negative charge
Ee       =>  nardi-pink         #negative charge
Ff       =>  nardi-dull-blue    #large hydrophobic
Gg       =>  nardi-orange       #glycine
Hh       =>  nardi-dull-blue    #large hydrophobic
Ii       =>  nardi-dull-blue    #hydrophobic
Kk       =>  nardi-red          #positive charge
Ll       =>  nardi-dull-blue    #hydrophobic
Mm       =>  nardi-dull-blue    #hydrophobic
Nn       =>  nardi-green        #polar
Pp       =>  nardi-yellow       #proline
Qq       =>  nardi-green        #polar
Rr       =>  nardi-red          #positive charge
Ss       =>  nardi-green        #small alcohol
Tt       =>  nardi-green        #small alcohol
Vv       =>  nardi-dull-blue    #hydrophobic
Ww       =>  nardi-dull-blue    #large hydrophobic
Yy       =>  nardi-cyan         #large hydrophobic
Zz       =>  nardi-pink         #E or Q
Xx       =>  nardi-dark-gray    #unknown
?        =>  nardi-light-gray   #unknown
*        =>  black              #stop

[CLUSTAL]
#Protein: highlight amino acid physicochemical properties
.        =>  clustal-dark-gray  #wildcard/mismatch
Aa       =>  clustal-blue       #hydrophobic
Bb       =>  clustal-white      #D or N
Cc       =>  clustal-pink       #hydrophobic
Dd       =>  clustal-magenta    #negative charge
Ee       =>  clustal-magenta    #negative charge
Ff       =>  clustal-blue       #large hydrophobic
Gg       =>  clustal-orange     #glycine
Hh       =>  clustal-cyan       #large hydrophobic
Ii       =>  clustal-blue       #hydrophobic
Kk       =>  clustal-red        #positive charge
Ll       =>  clustal-blue       #hydrophobic
Mm       =>  clustal-blue       #hydrophobic
Nn       =>  clustal-green      #polar
Pp       =>  clustal-yellow     #proline
Qq       =>  clustal-green      #polar
Rr       =>  clustal-red        #positive charge
Ss       =>  clustal-green      #small alcohol
Tt       =>  clustal-green      #small alcohol
Vv       =>  clustal-blue       #hydrophobic
Ww       =>  clustal-blue       #large hydrophobic
Yy       =>  clustal-cyan       #large hydrophobic
Zz       =>  clustal-white      #E or Q
Xx       =>  clustal-dark-gray  #unknown
?        =>  clustal-light-gray #unknown
*        =>  black              #stop

[CCLUSTAL]
#Protein consensus: highlight equivalence class
.        ->  nardi-dark-gray    #unconserved
+        ->  nardi-red          #positive charge
-        ->  nardi-pink         #negative charge
a        ->  nardi-dull-blue    #aromatic
c        ->  nardi-purple       #charged
h        ->  nardi-dull-blue    #hydrophobic
l        ->  nardi-dull-blue    #aliphatic
o        ->  nardi-green        #alcohol
p        ->  nardi-green        #polar
s        ->  nardi-blue         #small
t        ->  nardi-blue         #turnlike
u        ->  nardi-blue         #tiny

[KXLIN]
#Kuang Lin's colour neural net derived scheme
.        ->  dark-gray       #wildcard/mismatch
Aa       =>  lin-A
Rr       =>  lin-R
Nn       =>  lin-N
Dd       =>  lin-D
Cc       =>  lin-C
Qq       =>  lin-Q
Ee       =>  lin-E
Gg       =>  lin-G
Hh       =>  lin-H
Ii       =>  lin-I
Ll       =>  lin-L
Kk       =>  lin-K
Mm       =>  lin-M
Ff       =>  lin-F
Pp       =>  lin-P
Ss       =>  lin-S
Tt       =>  lin-T
Ww       =>  lin-W
Yy       =>  lin-Y
Vv       =>  lin-V
Bb       =>  light-gray      #D or N
Zz       =>  light-gray      #E or Q
Xx       ->  light-gray      #unknown
?        ->  light-gray      #unknown
*        =>  black           #stop

[DSSP]
#DSSP protein secondary structure assignments
.   ->  light-gray      #loop/irregular
Cc  ->  light-gray      #loop/irregular
Gg  =>  orange          #3-10 helix
Hh  =>  bright-red      #alpha helix
Ii  =>  brown           #pi helix
Ee  =>  bright-blue     #beta and beta bulge
Bb  =>  dull-blue       #residue in isolated beta bridge
Tt  =>  dark-green      #H-bonded turn
Ss  =>  dark-gray       #bend

[FIND]
#colour successive 'find' pattern blocks
Aa       =>  find-A
Bb       =>  find-B
Cc       =>  find-C
Dd       =>  find-D
Ee       =>  find-E
Ff       =>  find-F
Gg       =>  find-G
Hh       =>  find-H
Ii       =>  find-I
Jj       =>  find-J
Kk       =>  find-K
Ll       =>  find-L
Mm       =>  find-M
Nn       =>  find-N
Oo       =>  find-O
Pp       =>  find-P
Qq       =>  find-Q
Rr       =>  find-R
Ss       =>  find-S
Tt       =>  find-T

##########################################################################
