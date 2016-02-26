# Copyright (C) 1997-2006 Nigel P. Brown
# $Id: Display.pm,v 1.28 2015/06/14 17:09:03 npb Exp $

###########################################################################
=pod

=head1 Bio::MView::Display

Generate a multiple alignment figure for sequence subranges with respect to
a parent sequence. Supports special formatting per subrange including
HTML. Formats the alignment as records of required alignment width to fit
the page.

=cut
###########################################################################
package Bio::MView::Display;

use Bio::MView::Display::Any;

use strict;
use vars qw(%Template %Known_Types $Default_Stream
	    $Default_Columns $Default_Separator $Default_Gap $Default_Pad
	    $Default_Overlap $Default_Overlap_Color $Default_Ruler
	    $Default_HTML $Default_Bold $Default_Label0
	    $Default_Label1 $Default_Label2 $Default_Label3 $Default_Label4
	    $Default_Label5 $Default_Label6 $Intercolumn_Space);

$Default_Stream        = \*STDOUT;  #output stream
$Default_Columns       = 80;        #sequence display width
$Default_Separator     = "\n";      #output record separator
$Default_Gap           = '-';       #sequence gap character
$Default_Pad           = '';        #sequence trailing space character
$Default_Overlap       = '%';       #sequence overlap character
$Default_Overlap_Color = '#474747'; #sequence overlap color if HTML
$Default_Ruler         = 'abs';     #ruler is absolute from start of sequence
$Default_HTML          = 0;         #produce HTML?
$Default_Bold          = 0;         #embolden some text
$Default_Label0        = 1;         #show label0
$Default_Label1        = 1;         #show label1
$Default_Label2        = 1;         #show label2
$Default_Label3        = 1;         #show label3
$Default_Label4        = 1;         #show label4
$Default_Label5        = 1;         #show label5
$Default_Label6        = 1;         #show label6

$Intercolumn_Space     = ' ';       #space between output column types

%Template =
    (
     'string'    => undef, #parent sequence string
     'length'    => 0,     #parent sequence length
     'start'     => undef, #forward start position of parent sequence
     'stop'      => undef, #forward stop position of parent sequence
     'object'    => undef, #ordered list of display objects
     'posnwidth' => 0,     #string width of left/right sequence position
     'labwidth0' => 0,     #string width of first title block
     'labwidth1' => 0,     #string width of second title block
     'labwidth2' => 0,     #string width of third title block
     'labwidth3' => 0,     #string width of fourth title block
     'labwidth4' => 0,     #string width of fifth title block
     'labwidth5' => 0,     #string width of sixth title block
     'labwidth6' => 0,     #string width of seventh title block
    );

%Known_Types =
    (
     'Ruler'    => 1,
     'Header'   => 1,
     'Sequence' => 1,
     'Subrange' => 1,
    );

sub new {
    my $type = shift;
    die "$type::new() missing arguments\n"  unless @_;
    my $sequence = shift;
    my $self = { %Template };

    $self->{'string'}    = $sequence;
    $self->{'length'}    = $sequence->length;
    $self->{'start'}     = $sequence->lo;
    #$self->{'stop'}      = $sequence->hi;
    $self->{'stop'}      = $self->{'start'} + $self->{'length'} - 1;

    #my $tmp = abs($self->{'stop'}-$self->{'start'})+1;
    #warn "$self->{'start'} $self->{'stop'} :  $self->{length}  $tmp";

    $self->{'posnwidth'} = length(max($self->{'start'}, $self->{'stop'}));
    $self->{'object'}    = [];

    bless $self, $type;

    $self->append(@_)  if @_;

    $self;
}

#sub DESTROY { warn "DESTROY $_[0]\n" }

#force break of back references to allow garbage collection
sub free {
    local $_;
    #warn "FREE $_[0]\n";
    while ($_ = shift @{$_[0]->{'object'}}) { $_->free }
}

sub set_parameters {
    my $self = shift;
    my ($key, $val);

    while (@_) {
	($key, $val) = (shift, shift);
	if (exists $Template{$key}) {
	    $self->{$key} = $val;
	    next;
	}
	warn "Bio::MView::Manager: unknown parameter '$key'\n";
    }

    $self;
}

sub append {
    my $self = shift;
    #warn "${self}::append(@_)\n";
    my ($row, $type, $wide0, $wide1, $wide2, $wide3, $wide4, $wide5, $wide6);

    foreach $row (@_) {

	#check row 'type' exists and is valid
	if (exists $row->{'type'}) {
	    $type = ucfirst $row->{'type'};
	} else {
	    warn "$self append() missing type in '$row'\n";
	    next;
	}
	unless (exists $Known_Types{$type}) {
	    warn "$self append() unknown alignment type '$type'\n";
	    next;
	}

	#construct row object
	no strict 'refs';
	$type = "Bio::MView::Display::$type"->new($self, $row);
	use strict 'refs';

	push @{$self->{'object'}}, $type;

	$wide0 = length($type->label0);
	$wide1 = length($type->label1);
	$wide2 = length($type->label2);
	$wide3 = length($type->label3);
	$wide4 = length($type->label4);
	$wide5 = length($type->label5);
	$wide6 = length($type->label6);

	#warn $type->label0, ",", $type->label1, ",", $type->label2, "\n";
	
	$self->{'labwidth0'} = max($self->{'labwidth0'}, $wide0);
	$self->{'labwidth1'} = max($self->{'labwidth1'}, $wide1);
	$self->{'labwidth2'} = max($self->{'labwidth2'}, $wide2);
	$self->{'labwidth3'} = max($self->{'labwidth3'}, $wide3);
	$self->{'labwidth4'} = max($self->{'labwidth4'}, $wide4);
	$self->{'labwidth5'} = max($self->{'labwidth5'}, $wide5);
	$self->{'labwidth6'} = max($self->{'labwidth6'}, $wide6);
    }
    #Universal::vmstat("Display::append done");
    $self;
}

sub print {
    my $self = shift;
    local $_;
    foreach (sort keys %$self) {
	printf "%15s => %s\n", $_, $self->{$_};
    }
    map { $_->print } @{$self->{'object'}};
    $self;
}

sub length { $_[0]->{'length'} }

sub set_widths {
    my $self = shift;
    my ($key, $val);
    while (@_) {
	($key, $val) = (shift, shift);
	$self->{$key} = $val;
    }
    $self;
}

sub display {
    my $self = shift;
    my %par = @_;

    $par{'stream'} = $Default_Stream           unless exists $par{'stream'};
    $par{'col'}    = $Default_Columns  	       unless exists $par{'col'};
    $par{'gap'}    = $Default_Gap      	       unless exists $par{'gap'};
    $par{'pad'}    = $Default_Pad      	       unless exists $par{'pad'};
    $par{'rec'}    = $Default_Separator        unless exists $par{'rec'};
    $par{'ruler'}  = $Default_Ruler            unless exists $par{'ruler'};
    $par{'bold'}   = $Default_Bold             unless exists $par{'bold'};
    $par{'label0'} = $Default_Label0           unless exists $par{'label0'};
    $par{'label1'} = $Default_Label1           unless exists $par{'label1'};
    $par{'label2'} = $Default_Label2           unless exists $par{'label2'};
    $par{'label3'} = $Default_Label3           unless exists $par{'label3'};
    $par{'label4'} = $Default_Label4           unless exists $par{'label4'};
    $par{'label5'} = $Default_Label5           unless exists $par{'label5'};
    $par{'label6'} = $Default_Label6           unless exists $par{'label6'};

    $par{'posnwidth'} = $self->{'posnwidth'}   unless exists $par{'posnwidth'};

    $par{'labwidth0'} = $self->{'labwidth0'}   unless exists $par{'labwidth0'};
    $par{'labwidth1'} = $self->{'labwidth1'}   unless exists $par{'labwidth1'};
    $par{'labwidth2'} = $self->{'labwidth2'}   unless exists $par{'labwidth2'};
    $par{'labwidth3'} = $self->{'labwidth3'}   unless exists $par{'labwidth3'};
    $par{'labwidth4'} = $self->{'labwidth4'}   unless exists $par{'labwidth4'};
    $par{'labwidth5'} = $self->{'labwidth5'}   unless exists $par{'labwidth5'};
    $par{'labwidth6'} = $self->{'labwidth6'}   unless exists $par{'labwidth6'};

    $par{'html'} = 1    if $par{'bold'};
    if ($par{'html'}) {
	$par{'lap'}   = $Default_Overlap_Color unless exists $par{'lap'};
    } else {
	$par{'lap'}   = $Default_Overlap       unless exists $par{'lap'};
    }

    if ($par{'col'} < 1) {
	#full length display if meaningless column count
	$par{'col'} = $self->{'length'};
    }
    if (CORE::length $par{'gap'} != 1) {
	warn "${self}::display() gap must be single character '$par{'gap'}'\n";
	$par{'gap'} = $Default_Gap;
    }
    if (CORE::length $par{'pad'} > 1) {
	warn "${self}::display() pad must be null or single character '$par{'pad'}'\n";
	$par{'pad'} = $Default_Pad;
    }

    #map { warn "$_ => '$par{$_}'\n" } sort keys %par;

    my $str = $par{'stream'};

    my ($prefix, $suffix) = $par{'bold'} ? ('<STRONG>','</STRONG>') : ('','');

    print $str "<PRE>\n"    if $par{'html'};

    map { $_->reset } @{$self->{'object'}};

    my ($posnwidth, $o, $s, @left, @middle, @right, @tmp, $tmp) = (0);

    #need space for sequence numbers?
    foreach $o (@{$self->{'object'}}) {
	if ($o->{'number'}) {
	    $posnwidth = 1;
	    last;
	}
    }

    #iterate over display rows
LOOP:
    {
	while (1) {

	    last LOOP   unless scalar @{$self->{'object'}};

	    #Universal::vmstat("display LOOP");

	    #do record
	    foreach $o (@{$self->{'object'}}) {

		@left = @middle = @right = ();

		#do line
		$s = $o->next($par{'html'}, $par{'bold'}, $par{'col'},
			      $par{'gap'}, $par{'pad'}, $par{'lap'},
			      $par{'ruler'});

		last LOOP   unless $s;
		#warn "[ @$s ]\n";

		#### left ############
		@tmp = ();

		#first label
		if ($par{'label0'} and $par{'labwidth0'}) {
		    $tmp = $o->label0;
		    if ($tmp =~ /^\d+$/) {
			#numeric label - right justify
			push @tmp, sprintf("%$par{'labwidth0'}s", $tmp);
		    } else {
			#string label - left justify
			push @tmp, sprintf("%-$par{'labwidth0'}s", $tmp);
		    }
		    push @tmp, $Intercolumn_Space;
		}

		#second label
		if ($par{'label1'} and $par{'labwidth1'}) {
		    if ($par{'html'} and $o->{'url'}) {
			push @tmp,
			"<A HREF=\"$o->{'url'}\">", $o->label1, "</A>";
			push @tmp,
			" " x ($par{'labwidth1'}-CORE::length $o->label1);
		    } else {
			push @tmp,
			sprintf("%-$par{'labwidth1'}s", $o->label1);
		    }
		    push @tmp, $Intercolumn_Space;
		}

		#third label
		if ($par{'label2'} and $par{'labwidth2'}) {
		    push @tmp, sprintf("%-$par{'labwidth2'}s", $o->label2);
		    push @tmp, $Intercolumn_Space;
		}

		#fourth label
		if ($par{'label3'} and $par{'labwidth3'}) {
		    push @tmp, sprintf("%$par{'labwidth3'}s", $o->label3);
		    push @tmp, $Intercolumn_Space;
		}

		#fifth label
		if ($par{'label4'} and $par{'labwidth4'}) {
		    push @tmp, sprintf("%$par{'labwidth4'}s", $o->label4);
		    push @tmp, $Intercolumn_Space;
		}

		#sixth label
		if ($par{'label5'} and $par{'labwidth5'}) {
		    push @tmp, sprintf("%$par{'labwidth5'}s", $o->label5);
		    push @tmp, $Intercolumn_Space;
		}

		#seventh label
		if ($par{'label6'} and $par{'labwidth6'}) {
		    push @tmp, sprintf("%$par{'labwidth6'}s", $o->label6);
		    push @tmp, $Intercolumn_Space;
		}

		#left sequence position
		if ($posnwidth) {
		    if ($o->{'number'}) {
			push @tmp, sprintf("$prefix%$par{'posnwidth'}d$suffix", $s->[0]);
		    } else {
			push @tmp, sprintf("$prefix%$par{'posnwidth'}s$suffix", ' ');
		    }
		}

		push @tmp, $Intercolumn_Space;

		push @left, @tmp;

		#### middle ############
		@tmp = ();

		#sequence string
		push @tmp, $s->[2];
		push @tmp, $Intercolumn_Space;

		push @middle, @tmp;

		#### right ############
		@tmp = ();

		#right position
		if ($posnwidth) {
  		    if ($o->{'number'}) {
			push @tmp,
			sprintf("$prefix%-$par{'posnwidth'}d$suffix",
				$s->[1]);
		    } else {
			push @tmp,
			sprintf("$prefix%-$par{'posnwidth'}s$suffix", '');
		    }
		}

		push @tmp, "\n";

		push @right, @tmp;

		#### output row ############
		print $str @left, @middle, @right;

	    }#foreach

	    #record separator
	    print $str $par{'rec'}    unless $self->{'object'}->[0]->done;

	    #Universal::vmstat("display loop done");
 	}
    }
    print $str "</PRE>\n"    if $par{'html'};

    $self;
}

sub displaytext {
    return $_[0]    if defined $_[0] and $_[0] ne '';
    return '';
}

sub max { $_[0] > $_[1] ? $_[0] : $_[1] }


###########################################################################
1;

=pod

=head1 Synopsis

  my $sequence    #Bio::MView::Sequence instance
  my $from        #integer: start position of sequence for labelling
  my @list        #list of display rows

  #build an instance of an alignment object (mandatory arguments)

  my $a = new Bio::MView::Display($sequence);

  #build an instance of an alignment object (with extra parameters)

  my $b = new Bio::MView::Display($sequence, $from, @list);

  #add some display items to the $a object

  $a->append(@list);

  #explicitly add some (more?) display items to the $a object

  $a->append
      (

       { 'type'  => 'ruler' },

       { 'type'  => 'sequence',
	 'label1'=> 'first',
	 'url'   => 'http://fred.bloggs2',
	 'label2'=> 1.111,
	 'number'=> 1,
	 'range' => [ 1,  10, 'color'=>'red' ],
       },

       { 'type'  => 'subrange',
	 'label1'=> 'second',
	 'label2'=> 2.2,
	 'url'   => 'http://fred.bloggs',
	 'range' => [
		       1,  10, 'case'=>'uc','prefix'=>'<I>','suffix'=>'</I>',
		      15,  20, 'sym'=>'x','url'=>"URL1",
		      55, 104, 'sym'=>'*','url'=>"URL2",
		     100, 112, 'sym'=>'@','url'=>"URL3",
		    ],
       },

       { 'type'  => 'subrange',
	 'label1'=> 'third',
	 'label2'=> 3.14,
	 'range' => [
		      20, 100, 'color'=>'green' ],
		     500, 670, 'sym'=>'*','url'=>"URL4",
		     470, 510, 'color'=>'salmon',
		     800,1000, 'sym'=>'@',
		    ],
       },

       { 'type'  => 'subrange',
	 'label1'=> 'more data',
	 'label2'=> 4,
	 'range' => [
		      12, 314, 'sym'=>'h','color'=>'red',
		     200, 250, 'sym'=>'h','color'=>'red',
		     260, 260, 'sym='>'#',  #these two lines are equivalent
		          260, 'sym='>'#',  #this saves memory
		     328, 457, 'color'=>'purple',
		    ],
       },

      );

  #print the $a object using default formatting

  print $a->display;

  #print the $a object with some formatting overridden

  print $a->display('gap'=>'.', 'col'=>60, 'pad'=>' ', 'rec'=>"\n\n\n");


=head1 Defining a display list

A display list comprises a list of display items ordered as per the output
desired by the caller. Each display item may contain individual formatting
instructions for each subrange within. Currently, there are three types of
display item given by the 'type' attribute (see below).


 Attribute     Description
 -----------------------------------------------------------------------

 type          Mandatory. 'sequence' or 'subrange' or 'ruler', specifying
               the display of the parent sequence, some subranges thereof,
               or a simple ruler. Rulers count off 5, 10, 100's.

 label0        Optional. First label string of display item.

 label1        Optional. Second label string of display item.

 label2        Optional. Third label string of display item.

 label3        Optional. Fourth label string of display item.

 label4        Optional. Fifth label string of display item.

 label5        Optional. Sixth label string of display item.

 label6        Optional. Seventh label string of display item.

 url           Optional. URL to associate with 'label1'.

 number        Optional. If 1, attach sequence numbers to this item.

 sequence      Optional. Use supplied sequence instead of parent sequence
               as source of displayed characters. Must be same length as
               parent sequence.

 range         Optional. Any subranges of the parent sequence to be
               displayed for this display item, using optional formatting
               control attributes per subrange.

 paint         If paint=1 (default=0) then overlapping ranges for this row
               will use the colouring or symbol or URL of the last range
               definition in the region of the overlap as passed by the
               caller. The default behaviour is not to paint and to use the
               overlap symbol/colour (see below) or the last overlapping
               subrange by sequence order.

 prefix        if html, prefix whole row with this string.
 suffix        if html, suffix whole row with this string.

 pad           Pad character to fill leading/trailing empty space: overrides
               the display() pad setting for just this sequence.

The 'range' attribute is an anonymous array whose elements are also
anonymous arrays giving range and formatting information. Ranges can be
given in any order.  Each range consists minimally of a pair of ordered
integers indexing the parent sequence from 1, even if the actual labelling
for the parent sequence was specified differently.  An optional third
element after the range integer pair is an anonymous hash specifying
formatting details described below.



=head1 Formatting control: in display() call.

Display format control affecting the whole output is limited to:

 Attribute  Default      Behaviour
 -----------------------------------------------------------------------

 html       1            Controls HTML formatting {0=off,1=on}.

 bold       0            Embolden alignment

 col        80           Alignment width as number of columns. Setting
                         col < 1 causes the whole alignment to be
                         formatted as one long scrollable record.

 gap        '-'          Gap character to use.

 lap        '%'          Set the overlap character.
 lap        '#474747'    Set the overlap colour if html=1 and lap is more
                         than one character in length. Or set the overlap
                         character if html=1 and lap is just a single
                         character.

 pad        ''           Pad character to fill leading/trailing empty space.

 rec        "\n"         Record separator.

 ruler      'abs'        Controls ruler numbering: 'abs'= absolute from the
                         extreme left edge of the alignment, even if there
                         are marginal gaps. No other modes supported yet.

Overlaps can occur between subranges. If producing HTML, the display()
method will output the parent sequence in the overlap painted in the 'lap'
colour (defaulting to '#474747' [gray28] as above) or using the lap character
(defaulting to '%' as above). In non-HTML mode only the 'lap' character
will be printed. Clearly, if subranges overlap completely, all information
(colours, symbols, URLs) will be lost by this strategy. Alternatively, you
can simply ignore the cause of the overlap and continue with the new symbol
or colour or URL definitions by setting paint=1 for the whole row.  This is
useful for building up a complex image, for example, to define a long
subrange to be grey, then overlay some other colours at selected positions.



=head1 Formatting control: in subrange definitions.

Subrange format controls are:

 Attribute          Behaviour
 -------------------------------------------------------------------

 url                Embed range within this URL.

 sym                Use this symbol. Defaults to the parent sequence.

 case               Change range to case {uc|lc}. Ignored if 'sym' set.

 color              Set SPAN with this colour.

 class              Use SPAN.CLASS for this symbol.

 prefix		    Prefix with this arbitrary piece of HTML.

 suffix             Suffix with this arbitrary piece of HTML.



=head1 Example output

 <PRE>
		     [   .    :    .    :    .    :    .    :    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111   1 <SPAN style="color:red">GGAAATAAAA</SPAN>ATATAAAACCCGGTCGCAAGCGTGTCGCTGTAGATGAGCTGAAACGCAATTTTTCAGTGACTTTTTATCT  80
 <A HREF="http://fred.bloggs">second</A>      2.2     <I>GGAAATAAAA</I>----<A HREF="URL1">xxxxxx</A>----------------------------------<A HREF="URL2">**************************</A>
 third      3.14     -------------------<SPAN style="color:green">y5923y98y23v9t2 t oi2ltl3ntli2tit2309tu2905t2390tuvn 2u m t[3</SPAN>
 more data     4     -----------<SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN>

			 .    :    .    |    .    :    .    :    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111  81 CTCTAAAGACGAGCGTAGAAAGCGTCAATTCTTTTGTCAAACGCCACATTTTAAAAACAATCATTTACAAAAAAGGCACT 160
 <A HREF="http://fred.bloggs">second</A>      2.2     <A HREF="URL2">*******************</A><SPAN style="color:#474747">AAGCG</SPAN><A HREF="URL3">@@@@@@@@</A>------------------------------------------------
 third      3.14     <SPAN style="color:green"> mv[ n[ 2  2   205g5</SPAN>------------------------------------------------------------
 more data     4     <SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN>

			 .    :    .    :    .    :    .    |    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111 161 AACXAAGATTCTTCTATCAATTGTGATTCTTCTAGTAGGCTTTAAGCCTACTTTAAGGTATCAGCTCTTAAGAGCAGCAC 240
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     --------------------------------------------------------------------------------
 more data     4     <SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN><SPAN style="color:#474747">CTTTAAGCCTACTTTAAGGTATCAGCTCTTAAGAGCAGCAC</SPAN>

			 .    :    .    :    .    :    .    :    .    :    .    |    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111 241 CCGCTAGAGCATTTAGTTTTAXTCCCTACCGATGAGCGATCCATCAAGCAAAACCAAACCCGGCCACAAACAAGCCTCCA 320
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     --------------------------------------------------------------------------------
 more data     4     <SPAN style="color:#474747">CCGCTAGAGC</SPAN><SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN>------

			 .    :    .    :    .    :    .    :    .    :    .    :    .    :    .    |
 <A HREF="http://fred.bloggs2">first</A>     1.111 321 GCAAATAACGCXAAGTCTTTAACCACTAGCCAAAATAGCTATATGCATCAAATAGCCGCCTAAATTTTGAAGTTTTGTAA 400
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     --------------------------------------------------------------------------------
 more data     4     -------<SPAN style="color:purple">ACGCXAAGTCTTTAACCACTAGCCAAAATAGCTATATGCATCAAATAGCCGCCTAAATTTTGAAGTTTTGTAA</SPAN>

			 .    :    .    :    .    :    .    :    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111 401 TCACTXCAAGCAATGATTTTAACGCTTGCATGGTTTTATTCCTTTAGGATATTTTTTATATTGCTAAAAACAGCCTAGAA 480
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     ---------------------------------------------------------------------<SPAN style="color:salmon">GFKGFKFHFHF</SPAN>
 more data     4     <SPAN style="color:purple">TCACTXCAAGCAATGATTTTAACGCTTGCATGGTTTTATTCCTTTAGGATATTTTTT</SPAN>-----------------------

			 .    :    .    |    .    :    .    :  ]
 <A HREF="http://fred.bloggs2">first</A>     1.111 481 ATGATAGCGCATTATGTTTAATTXCAACAAAATCATTACTATT 523
 <A HREF="http://fred.bloggs">second</A>      2.2     -------------------------------------------
 third      3.14     <SPAN style="color:salmon">YTRU7r7r76576565576</SPAN><SPAN style="color:#474747">77%%AGCGCAT</SPAN><A HREF="URL4">*************</A>
 more data     4     -------------------------------------------

 </PRE>


=head1 Example output formatted

 <PRE>
		     [   .    :    .    :    .    :    .    :    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111   1 <SPAN style="color:red">GGAAATAAAA</SPAN>ATATAAAACCCGGTCGCAAGCGTGTCGCTGTAGATGAGCTGAAACGCAATTTTTCAGTGACTTTTTATCT  80
 <A HREF="http://fred.bloggs">second</A>      2.2     <I>GGAAATAAAA</I>----<A HREF="URL1">xxxxxx</A>----------------------------------<A HREF="URL2">**************************</A>
 third      3.14     -------------------<SPAN style="color:green">y5923y98y23v9t2 t oi2ltl3ntli2tit2309tu2905t2390tuvn 2u m t[3</SPAN>
 more data     4     -----------<SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN>

			 .    :    .    |    .    :    .    :    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111  81 CTCTAAAGACGAGCGTAGAAAGCGTCAATTCTTTTGTCAAACGCCACATTTTAAAAACAATCATTTACAAAAAAGGCACT 160
 <A HREF="http://fred.bloggs">second</A>      2.2     <A HREF="URL2">*******************</A><SPAN style="color:#474747">AAGCG</SPAN><A HREF="URL3">@@@@@@@@</A>------------------------------------------------
 third      3.14     <SPAN style="color:green"> mv[ n[ 2  2   205g5</SPAN>------------------------------------------------------------
 more data     4     <SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN>

			 .    :    .    :    .    :    .    |    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111 161 AACXAAGATTCTTCTATCAATTGTGATTCTTCTAGTAGGCTTTAAGCCTACTTTAAGGTATCAGCTCTTAAGAGCAGCAC 240
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     --------------------------------------------------------------------------------
 more data     4     <SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN><SPAN style="color:#474747">CTTTAAGCCTACTTTAAGGTATCAGCTCTTAAGAGCAGCAC</SPAN>

			 .    :    .    :    .    :    .    :    .    :    .    |    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111 241 CCGCTAGAGCATTTAGTTTTAXTCCCTACCGATGAGCGATCCATCAAGCAAAACCAAACCCGGCCACAAACAAGCCTCCA 320
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     --------------------------------------------------------------------------------
 more data     4     <SPAN style="color:#474747">CCGCTAGAGC</SPAN><SPAN style="color:red">hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh</SPAN>------

			 .    :    .    :    .    :    .    :    .    :    .    :    .    :    .    |
 <A HREF="http://fred.bloggs2">first</A>     1.111 321 GCAAATAACGCXAAGTCTTTAACCACTAGCCAAAATAGCTATATGCATCAAATAGCCGCCTAAATTTTGAAGTTTTGTAA 400
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     --------------------------------------------------------------------------------
 more data     4     -------<SPAN style="color:purple">ACGCXAAGTCTTTAACCACTAGCCAAAATAGCTATATGCATCAAATAGCCGCCTAAATTTTGAAGTTTTGTAA</SPAN>

			 .    :    .    :    .    :    .    :    .    :    .    :    .    :    .    :
 <A HREF="http://fred.bloggs2">first</A>     1.111 401 TCACTXCAAGCAATGATTTTAACGCTTGCATGGTTTTATTCCTTTAGGATATTTTTTATATTGCTAAAAACAGCCTAGAA 480
 <A HREF="http://fred.bloggs">second</A>      2.2     --------------------------------------------------------------------------------
 third      3.14     ---------------------------------------------------------------------<SPAN style="color:salmon">GFKGFKFHFHF</SPAN>
 more data     4     <SPAN style="color:purple">TCACTXCAAGCAATGATTTTAACGCTTGCATGGTTTTATTCCTTTAGGATATTTTTT</SPAN>-----------------------

			 .    :    .    |    .    :    .    :  ]
 <A HREF="http://fred.bloggs2">first</A>     1.111 481 ATGATAGCGCATTATGTTTAATTXCAACAAAATCATTACTATT 523
 <A HREF="http://fred.bloggs">second</A>      2.2     -------------------------------------------
 third      3.14     <SPAN style="color:salmon">YTRU7r7r76576565576</SPAN><SPAN style="color:#474747">77%%AGCGCAT</SPAN><A HREF="URL4">*************</A>
 more data     4     -------------------------------------------

 </PRE>


=head1 Copyright (C) 1997-2006 Nigel P. Brown.

$Id: Display.pm,v 1.28 2015/06/14 17:09:03 npb Exp $


=cut
