#!/usr/bin/env perl
# Copyright (C) 2015 Nigel P. Brown

use Test::More;
use lib 'lib';
use strict;

print "test__loading\n";
use_ok('Bio::MView::Sequence');

sub main {
    my @args = @_;

    test__create_empty();
    test__forward_insert();
    test__reverse_insert();
    test__special_characters();
    test__toy_cases();

    done_testing;
    return 0;
}

exit main(@ARGV);

###########################################################################
TODO: { local $TODO = "not implemented"; };

#used repeatedly by other tests
sub test_state {
    my ($S, $h) = @_;
    if ($h->{reverse}) {
	isa_ok($S, 'Bio::MView::Reverse_Sequence');
	is($S->is_reversed, $h->{reverse}, "is_reverse() true");
    } else {
	isa_ok($S, 'Bio::MView::Sequence');
	is($S->is_reversed, $h->{reverse}, "is_reverse() false");
    }
    is($S->frameshifts, $h->{fs},	  "frameshifts() is $h->{fs}");
    is($S->lo, $h->{lo},		  "lo() is $h->{lo}");
    is($S->hi, $h->{hi},		  "hi() is $h->{hi}");
    is($S->reflo, $h->{reflo},      	  "reflo() is $h->{reflo}");
    is($S->refhi, $h->{refhi},       	  "refhi() is $h->{refhi}");
    is($S->leader, $h->{leader},          "leader() is $h->{leader}");
    is($S->trailer, $h->{trailer}, 	  "trailer() is $h->{trailer}");
    is($S->trulen, $h->{trulen}, 	  "trulen() is $h->{trulen}");
    is($S->lablen, $h->{lablen},          "lablen() is $h->{lablen}");
    #is($S->length, $h->{length}, 	  "length() is $h->{length}");
    is($S->seqlen, $h->{seqlen},          "seqlen() is $h->{seqlen}");
    is($S->string, $h->{string},	  "string() is '$h->{string}'");
    is($S->sequence, $h->{sequence},      "sequence() is '$h->{sequence}'");
}

sub test_labels {
    my ($S, $h) = @_;
    is($S->fromlabel1, $h->{fromlabel1}, "fromlabel1() is $h->{fromlabel1}");
    is($S->tolabel1, $h->{tolabel1},	 "tolabel1() is $h->{tolabel1}");
    is($S->fromlabel2, $h->{fromlabel2}, "fromlabel2() is $h->{fromlabel2}");
    is($S->tolabel2, $h->{tolabel2},     "tolabel2() is $h->{tolabel2}");
}

sub test_special_chars {
    my ($S, $h) = @_;
    is($S->get_pad, $h->{pad},	"get_pad() is '$h->{pad}'");
    is($S->get_gap, $h->{gap},	"get_gap() is '$h->{gap}'");
    is($S->get_spc, $h->{spc},	"get_spc() is '$h->{spc}'");
    is($S->get_fs1, $h->{fs1},	"get_fs1() is '$h->{fs1}'");
    is($S->get_fs2, $h->{fs2},	"get_fs2() is '$h->{fs2}'");
}

###########################################################################
sub test__create_empty {
    print "\ntest__create_empty\n";
    my $h = {
	'lo'       => 0,     'hi'          => 0,
	'reflo'    => 0,     'refhi'       => 0,
	'leader'   => 0,     'trailer'     => 0,
	'string'   => '',    'length'      => 0,
	'trulen'   => 0,     'lablen'      => 0,
	'sequence' => '',    'seqlen'      => 0,
	'reverse'  => undef, 'fs' => 0,
    };
    my $S;

    print "\ntest__create_empty / forward\n";
    $S = new Bio::MView::Sequence;
    ok(defined $S, "new Sequence() defined");
    $h->{reverse} = 0; test_state($S, $h);
    
    print "\ntest__create_empty / reverse\n";
    $S = new Bio::MView::Reverse_Sequence;
    ok(defined $S, "new Reverse_Sequence() defined");
    $h->{reverse} = 1; test_state($S, $h);

    print "\ntest__create_empty / default special characters\n";
    my $c = {
	'pad' => qw(.), 'gap' => qw(-), 'spc' => ' ',
	'fs1' => qw(/), 'fs2' => qw(\\),
    };
    test_special_chars($S, $c);
}

sub test__forward_insert {
    print "\ntest__forward_insert: B AB ABC X-ABC X-ABC-Y\n";
    my ($S, $s, $f, $h);

    print "\ntest__forward_insert / insert empty\n";
    $S = new Bio::MView::Sequence();
    $h = {
	'lo'       => 0,    	'hi'             => 0,
	'reflo'    => 0,    	'refhi'          => 0,
	'leader'   => 0,    	'trailer'        => 0,
	'string'   => '',    	'length'         => 0,
	'trulen'   => 0,        'lablen'         => 0,
	'sequence' => '',    	'seqlen'         => 0,
	'reverse'  => 0,    	'fs' => 0,
    };
    ok($S->insert(), "insert(empty)"); test_state($S, $h);
    #end

    print "\ntest__forward_insert / insert character at [20]\n";
    $S = new Bio::MView::Sequence();
    $s = 'B'; $f = [\$s, 20, 20]; $h = {
        'lo'       => 20,   	'hi'             => 20,
        'reflo'    => 20,   	'refhi'          => 20,
        'leader'   => 0,    	'trailer'        => 0,
        'string'   => 'B',  	'length'         => 1,
	'trulen'   => 1,        'lablen'         => 1,
        'sequence' => 'B',  	'seqlen'         => 1,
        'reverse'  => 0,    	'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__forward_insert / overwrite character at [20]\n";
    $S->insert($f); test_state($S, $h);

    print "\ntest__forward_insert / prepend character at [19]\n";
    $s = 'A'; $f = [\$s, 19, 19]; $h = {
        'lo'       => 19,       'hi'             => 20,
        'reflo'    => 19,       'refhi'          => 20,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => 'AB',     'length'         => 2,
	'trulen'   => 2,        'lablen'         => 2,
        'sequence' => 'AB',     'seqlen'         => 2,
        'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__forward_insert / insert character at [21]\n";
    $s = 'C'; $f = [\$s, 21, 21]; $h = {
        'lo'       => 19,       'hi'             => 21,
        'reflo'    => 19,       'refhi'          => 21,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => 'ABC',    'length'         => 3,
	'trulen'   => 3,        'lablen'         => 3,
        'sequence' => 'ABC',    'seqlen'         => 3,
        'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__forward_insert / prepend character at [17] with a gap\n";
    $s = 'X'; $f = [\$s, 17, 17]; $h = {
        'lo'       => 17,       'hi'             => 21,
        'reflo'    => 17,       'refhi'          => 21,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => 'X-ABC',  'length'         => 5,
	'trulen'   => 5,        'lablen'         => 5,
        'sequence' => 'XABC',   'seqlen'         => 4,
        'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__forward_insert / apppend character at [23] with a gap\n";
    $s = 'Y'; $f = [\$s, 23, 23]; $h = {
        'lo'       => 17,        'hi'             => 23,
        'reflo'    => 17,        'refhi'          => 23,
        'leader'   => 0,         'trailer'        => 0,
        'string'   => 'X-ABC-Y', 'length'         => 7,
	'trulen'   => 7,         'lablen'         => 7,
        'sequence' => 'XABCY',   'seqlen'         => 5,
        'reverse'  => 0,         'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end
}

sub test__reverse_insert {
    print "\ntest__reverse_insert: B AB ABC X-ABC X-ABC-Y\n";
    my ($S, $s, $f, $h);

    print "\ntest__reverse_insert / insert empty\n";
    $S = new Bio::MView::Reverse_Sequence();
    $h = {
        'lo'       => 0,        'hi'             => 0,
        'reflo'    => 0,        'refhi'          => 0,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => '',       'length'         => 0,
        'trulen'   => 0,        'lablen'         => 0,
        'sequence' => '',       'seqlen'         => 0,
        'reverse'  => 1,        'fs' => 0,
    };
    ok($S->insert(), "insert(empty)"); test_state($S, $h);
    #end

    print "\ntest__reverse_insert / insert character at [20]\n";
    $S = new Bio::MView::Reverse_Sequence();
    $s = 'B'; $f = [\$s, 20, 20]; $h = {
        'lo'       => 20,       'hi'             => 20,
        'reflo'    => 20,       'refhi'          => 20,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => 'B',      'length'         => 1,
        'trulen'   => 1,        'lablen'         => 1,
        'sequence' => 'B',      'seqlen'         => 1,
        'reverse'  => 1,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__reverse_insert / overwrite character at [20]\n";
    $S->insert($f); test_state($S, $h);

    print "\ntest__reverse_insert / prepend character at [19]\n";
    $s = 'A'; $f = [\$s, 19, 19]; $h = {
        'lo'       => 19,       'hi'             => 20,
        'reflo'    => 20,       'refhi'          => 19,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => 'BA',     'length'         => 2,
        'trulen'   => 2,        'lablen'         => 2,
        'sequence' => 'BA',     'seqlen'         => 2,
        'reverse'  => 1,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__reverse_insert / insert character at [21]\n";
    $s = 'C'; $f = [\$s, 21, 21]; $h = {
        'lo'       => 19,       'hi'             => 21,
        'reflo'    => 21,       'refhi'          => 19,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => 'CBA',    'length'         => 3,
        'trulen'   => 3,        'lablen'         => 3,
        'sequence' => 'CBA',    'seqlen'         => 3,
        'reverse'  => 1,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__reverse_insert / prepend character at [17] with a gap\n";
    $s = 'X'; $f = [\$s, 17, 17]; $h = {
        'lo'       => 17,       'hi'             => 21,
        'reflo'    => 21,       'refhi'          => 17,
        'leader'   => 0,        'trailer'        => 0,
        'string'   => 'CBA-X',  'length'         => 5,
        'trulen'   => 5,        'lablen'         => 5,
        'sequence' => 'CBAX',   'seqlen'         => 4,
        'reverse'  => 1,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__reverse_insert / apppend character at [23] with a gap\n";
    $s = 'Y'; $f = [\$s, 23, 23]; $h = {
        'lo'       => 17,        'hi'             => 23,
        'reflo'    => 23,        'refhi'          => 17,
        'leader'   => 0,         'trailer'        => 0,
        'string'   => 'Y-CBA-X', 'length'         => 7,
        'trulen'   => 7,         'lablen'         => 7,
        'sequence' => 'YCBAX',   'seqlen'         => 5,
        'reverse'  => 1,         'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end
}

sub test__special_characters {
    print "\ntest__special_characters\n";
    my ($S, $s, $f, $h);

    print "\ntest__special_characters / pad\n";
    $S = new Bio::MView::Sequence();
    $S->set_pad('~');
    $s = '.P..'; $f = [\$s, 20, 23]; $h = {
	'lo'       => 20,       'hi'             => 23,
	'reflo'    => 21,       'refhi'          => 21,
	'leader'   => 1,        'trailer'        => 2,
	'string'   => '~P~~',   'length'         => 4,
	'trulen'   => 4,        'lablen'         => 4,
	'sequence' => 'P',      'seqlen'         => 1,
	'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end

    print "\ntest__special_characters / gap\n";
    $S = new Bio::MView::Sequence();
    $S->set_gap('?');
    $s = 'A-B'; $f = [\$s, 20, 22]; $h = {
	'lo'       => 20,       'hi'             => 22,
	'reflo'    => 20,       'refhi'          => 22,
	'leader'   => 0,        'trailer'        => 0,
	'string'   => 'A?B',    'length'         => 3,
	'trulen'   => 3,        'lablen'         => 3,
	'sequence' => 'AB',     'seqlen'         => 2,
	'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end

    print "\ntest__special_characters / spc\n";
    $S = new Bio::MView::Sequence();
    $S->set_gap('?');
    $s = ' A B '; $f = [\$s, 20, 24]; $h = {
	'lo'       => 20,       'hi'             => 24,
	'reflo'    => 20,       'refhi'          => 24,
	'leader'   => 0,        'trailer'        => 0,
	'string'   => ' A B ',  'length'         => 5,
	'trulen'   => 5,        'lablen'         => 5,
	'sequence' => 'AB',     'seqlen'         => 2,
	'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end

    print "\ntest__special_characters / fs1\n";
    $S = new Bio::MView::Sequence();
    $s = 'PES/GRP'; $f = [\$s, 20, 27]; $h = {
        'lo'       => 20,        'hi'             => 27,
        'reflo'    => 20,        'refhi'          => 27,
        'leader'   => 0,         'trailer'        => 0,
        'string'   => 'PES/GRP', 'length'         => 'FIXME',
	'trulen'   => 7,         'lablen'         => 8, #WHY
        'sequence' => 'PESGRP',  'seqlen'         => 6,
        'reverse'  => 0,         'fs' => 1,
    };
    $S->insert($f); test_state($S, $h);
    #end

    print "\ntest__special_characters / fs2\n";
    $S = new Bio::MView::Sequence();
    $s = 'PES\GRP'; $f = [\$s, 20, 26]; $h = {
        'lo'       => 20,        'hi'             => 26,
        'reflo'    => 20,        'refhi'          => 26,
        'leader'   => 0,         'trailer'        => 0,
        'string'   => 'PES\GRP', 'length'         =>'FIXME',
	'trulen'   => 7,         'lablen'         => 7,
        'sequence' => 'PESGRP',  'seqlen'         => 6,
        'reverse'  => 0,         'fs' => 1,
    };
    $S->insert($f); test_state($S, $h);
    #end

    print "\ntest__special_characters / fastx_3.4t23.dat (simple)\n";
    my $frg = [1, 115, 'TDQLEDEKSALQTEIANLLKEKEKLEFILAAHRPACKIPDDLGFPEEMSVASLDLTGGLPEVATPESEEAFTLPLLNDPEPKPSVEPVKSISSMELKTEPFDDFLFPASSRPSG'];
    my $ungapped = $frg->[2]; $ungapped =~ s/\///g;
    $S = new Bio::MView::Sequence();
    $f = [\$frg->[2], $frg->[0], $frg->[1]]; $h = {
        'lo'       => 1,         'hi'             => 115,
        'reflo'    => 1,         'refhi'          => 115,
        'leader'   => 0,         'trailer'        => 0,
        'string'   => $frg->[2], 'length'         => 'FIXME',
	'trulen'   => 114,       'lablen'         => 115, #WHY
        'sequence' => $ungapped, 'seqlen'         => length($ungapped),
        'reverse'  => 0,         'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end

    print "\ntest__special_characters / fastx_3.4t23.dat (frameshift)\n";
    my $frg = [18, 194, 'PAEGEGKTRVHPGSSPTCLQDP**PGLPRRDVCGFP*SDWGPARGCHPGV*GGLHPASPQ*P*AQ/SPQWNLSRASAAWS*RPSPLMTSCSQHHPGPVALRQPAPCQTWTYLGPSMQQTGSLCTVAPWGWGPWPQSWSPCALRWSPVLPAALLTRLPSSSPTPRLTPSPAVQLPTA'];
    my $ungapped = $frg->[2]; $ungapped =~ s/\///g;
    $S = new Bio::MView::Sequence();
    $f = [\$frg->[2], $frg->[0], $frg->[1]]; $h = {
        'lo'       => 18,        'hi'             => 194,
        'reflo'    => 18,        'refhi'          => 194,
        'leader'   => 0,         'trailer'        => 0,
        'string'   => $frg->[2], 'length'         => 'FIXME',
	'trulen'   => 176,       'lablen'         => 177, #WHY
        'sequence' => $ungapped, 'seqlen'         => length($ungapped),
        'reverse'  => 0,         'fs' => 1,
    };
    $S->insert($f); test_state($S, $h);
    #end
}

sub test__toy_cases {
    print "\ntest__toy_cases\n";
    my ($S, $s, $f, $h);

    print "\ntest__case 1: forward query\n";
    $S = new Bio::MView::Sequence();
    $s = '23456'; $f = [\$s, 2, 6]; $h = {
	'lo'       => 2,        'hi'             => 6,
	'reflo'    => 2,        'refhi'          => 6,
	'leader'   => 0,        'trailer'        => 0,
	'string'   => '23456',  'length'         => 'FIXME',
	'trulen'   => 5,        'lablen'         => 5,
	'sequence' => '23456',  'seqlen'         => 5,
	'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__case 1: forward query + forward hit\n";
    $S = new Bio::MView::Sequence();
    $S->set_pad('-');
    $s = '-23--'; $f = [\$s, 2, 6]; $h = {
	'lo'       => 2,        'hi'             => 6,
	'reflo'    => 3,        'refhi'          => 4,
	'leader'   => 1,        'trailer'        => 2,
	'string'   => '-23--',  'length'         => 'FIXME',
	'trulen'   => 5,        'lablen'         => 5,
	'sequence' => '23',     'seqlen'         => 2,
	'reverse'  => 0,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end

    print "\ntest__case 1: reverse query\n";
    $S = new Bio::MView::Reverse_Sequence();
    $s = '65432'; $f = [\$s, 6, 2]; $h = {
	'lo'       => 2,        'hi'             => 6,
	'reflo'    => 6,        'refhi'          => 2,
	'leader'   => 0,        'trailer'        => 0,
	'string'   => '65432',  'length'         => 'FIXME',
	'trulen'   => 5,        'lablen'         => 5,
	'sequence' => '65432',  'seqlen'         => 5,
	'reverse'  => 1,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);

    print "\ntest__case 1: reverse query + forward hit\n";
    $S = new Bio::MView::Reverse_Sequence();
    $S->set_pad('-');
    $s = '-23--'; $f = [\$s, 6, 2]; $h = {
	'lo'       => 2,        'hi'             => 6,
	'reflo'    => 5,        'refhi'          => 4,
	'leader'   => 1,        'trailer'        => 2,
	'string'   => '-23--',  'length'         => 'FIXME',
	'trulen'   => 5,        'lablen'         => 5,
	'sequence' => '23',     'seqlen'         => 2,
	'reverse'  => 1,        'fs' => 0,
    };
    $S->insert($f); test_state($S, $h);
    #end
}




__DATA__;

#TODO

test raw, col  (easy - use the tests I just did on X-ABC-Y)

test longer segments and al classes of overlap

1. ---- ----

2. -----           3.       -----
      -----              -----

4. ------------    5.     ---- 
       ----           -------------

(6.   ------     ------   ------   and the 3 opposites)
(     ------     -----     -----                      )

test negative? positions: should/does it fail?

test overwrite really does keep original character(s)

test labels:
    labels1  for the reference sequence
    labels2  for the dependent sequence

test special characters insertion/removal, especially frameshifts

test is_X predicates

--

repeat everything for Reverse_Sequence

--

subclass a Forward_Sequence

--

test _substr

implement a true substr

--

(test findall) - should it be here?

