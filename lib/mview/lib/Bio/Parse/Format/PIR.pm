# Copyright (C) 1998-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

###########################################################################
package Bio::Parse::Format::PIR;

use Bio::Parse::Strings;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::Parse::Record);


#PIR record types
my $PIR_Null     = '^\s*$';#'
my $PIR_SEQ      = '^\s*>';
my $PIR_SEQend   = $PIR_SEQ;


#Consume one entry-worth of input on text stream associated with $file and
#return a new PIR instance.
sub get_entry {
    my $text = shift;
    my $line = '';
    my $data = 0;

    while ($text->getline(\$line)) {

        #start of entry
        if (!$data) {
            $text->start_count();
            $data = 1;
            next;
        }

    }
    return 0  unless $data;

    new Bio::Parse::Format::PIR(undef, $text, $text->get_start(), $text->get_stop());
}

#Parse one entry
sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    while (defined ($line = $scan->read_line)) {

        #SEQ lines
        if ($line =~ /$PIR_SEQ/o) {
            $scan->scan_until($PIR_SEQend);
            $self->push_record('SEQ',
                               $scan->get_block_start(),
                               $scan->get_block_stop(),
                );
            next;
        }

        #blank line or empty record: ignore
        next  if $line =~ /$PIR_Null/o;

        #default
        $self->warn("unknown field: $line");
    }
    $self;#->examine;
}


###########################################################################
package Bio::Parse::Format::PIR::SEQ;

use Bio::Parse::Strings qw(strip_leading_space strip_trailing_space);

use vars qw(@ISA);

@ISA = qw(Bio::Parse::Record);

sub new {
    my $self = new Bio::Parse::Record(@_);
    my $scan = new Bio::Parse::Scanner($self);
    my $line = '';

    $self->{'prefix'} = '';
    $self->{'id'}     = '';
    $self->{'desc'}   = '';
    $self->{'seq'}    = '';

    while (defined ($line = $scan->read_line(1))) {

        #read header line
        if ($line =~ /^\s*>\s*(..);(\S+)/o) {
            $self->test_args(\$line, $1, $2);
            (
             $self->{'prefix'},
             $self->{'id'},
            ) = ($1, $2);

            #force read of next line for description
            $self->{'desc'} = $scan->read_line(1);
            $self->{'desc'} = strip_leading_space($self->{'desc'});
            $self->{'desc'} = strip_trailing_space($self->{'desc'});

            next;
        }

        #read sequence lines upto asterisk, if present
        if ($line =~ /([^\*]+)/) {
            $self->{'seq'} .= $1;
            next;
        }

        #ignore lone asterisk
        last    if $line =~ /\*/;

        #default
        $self->warn("unknown field: $line");

    }
    #strip internal whitespace from sequence
    $self->{'seq'} =~ s/\s//g;

    $self;
}

sub dump_data {
    my ($self, $indent) = (@_, 0);
    my $x = ' ' x $indent;
    my $s = '';
    $s .= sprintf "$x%20s -> %s\n",   'prefix', $self->{'prefix'};
    $s .= sprintf "$x%20s -> %s\n",   'id',     $self->{'id'};
    $s .= sprintf "$x%20s -> '%s'\n", 'desc',   $self->{'desc'};
    $s .= sprintf "$x%20s -> %s\n",   'seq',    $self->{'seq'};
    return $s;
}


###########################################################################
1;
