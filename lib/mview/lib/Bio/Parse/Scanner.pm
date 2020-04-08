# Copyright (C) 1996-2019 Nigel P. Brown

# This file is part of MView.
# MView is released under license GPLv2, or any later version.

use strict;

###########################################################################
package Bio::Parse::Scanner;

use Bio::Util::Object;
use Bio::Parse::ReadText;

use vars qw(@ISA);

@ISA = qw(Bio::Util::Object);

sub new {
    my ($type, $entry, $indent,  $text) = (@_, 0, undef);

    my $self = {};
    bless $self, $type;

    if (defined $text and ref $text) {
        #use supplied text and positions (or defaults thereof)
        $self->{'text'}   = new Bio::Parse::ReadText($text);
        $self->{'offset'} = 0;
        $self->{'bytes'}  = length($$text);
    } else {
        #use entry object's text and positions
        $self->{'text'}   = $entry->get_text();
        $self->{'offset'} = $entry->get_offset();
        $self->{'bytes'}  = $entry->get_bytes();
    }

    $self->{'indent'}     = $indent;
    $self->{'extent'}     = $self->{'offset'} + $self->{'bytes'};

    $self->{'linestart'}  = $self->{'offset'};
    $self->{'blockstart'} = $self->{'offset'};
    $self->{'cursor'}     = $self->{'offset'};

    $self->{'line'}       = '';
    $self->{'block'}      = '';

    $self;
}

sub get_line_start  { $_[0]->{'linestart'}; }
sub get_line_stop   { $_[0]->{'cursor'}; }
sub get_line_bytes  { $_[0]->{'cursor'} - $_[0]->{'linestart'}; }
sub get_line        { $_[0]->{'line'}; }

sub get_block_start { $_[0]->{'blockstart'}; }
sub get_block_stop  { $_[0]->{'cursor'}; }
sub get_block_bytes { $_[0]->{'cursor'} - $_[0]->{'blockstart'}; }
sub get_block       { $_[0]->{'block'}; }

###########################################################################
# read operations
###########################################################################
# read next line; chomp line if optional argument is non-zero;
# return line read or undef if no bytes.
sub read_line {
    my ($self, $chomp) = (@_, 0);
    #warn "read_line(chomp=$chomp)\n";

    return undef  unless $self->_read_line;

    chomp $self->{'line'}  if $chomp;

    return $self->{'line'};
}

# Read $count lines; return block read.
sub read_lines {
    my ($self, $count) = (@_, 0);
    #warn "read_lines: $count\n";

    my $line  = \$self->{'line'};
    my $block = \$self->{'block'};

    $$block = $$line; $self->{'blockstart'} = $self->{'linestart'};

    while ($count > 0 and $self->_read_line) {
        $$block .= $$line;
        $count--;
    }

    return $$block;
}

# read all remaining lines until EOF; return block read.
sub read_remainder {
    my ($self, $count, $key) = (@_, 0);
    #warn "read_remainder\n";

    my $line  = \$self->{'line'};
    my $block = \$self->{'block'};

    $$block = $$line; $self->{'blockstart'} = $self->{'linestart'};

    while ($self->_read_line) {
        $$block .= $$line;
    }

    return $$block;
}

# read >= 1 record lines terminating on matching 'pattern';
# return block read.
sub read_until_inclusive {
    my ($self, $pattern) = @_;
    #warn "read_until_inclusive: /$pattern/\n";

    my $line  = \$self->{'line'};
    my $block = \$self->{'block'};

    $$block = $$line; $self->{'blockstart'} = $self->{'linestart'};

    while ($self->_read_line) {
        $$block .= $$line;
        last  if $$line =~ /$pattern/;
    }

    return $$block;
}

# read >= 1 record lines terminating on matching 'pattern';
# return block read.
sub read_while {
    my ($self, $pattern) = @_;
    #warn "scan_while: /$pattern/\n";

    my $line  = \$self->{'line'};
    my $block = \$self->{'block'};

    $$block = $$line; $self->{'blockstart'} = $self->{'linestart'};

    while ($self->_read_line) {
        if ($$line !~ /$pattern/) {  #backup
            $self->{'cursor'} = $self->{'linestart'};
            $$line = '';
            last;
        }
        $$block .= $$line;
    }

    return $$block;
}

###########################################################################
# scan operations
###########################################################################
# read >= 1 record lines terminating on matching 'pattern'.
sub scan_until {
    my ($self, $pattern) = @_;
    #warn "scan_until: /$pattern/\n";

    my $line  = \$self->{'line'};

    $self->{'blockstart'} = $self->{'linestart'};

    while ($self->_read_line) {
        if ($$line =~ /$pattern/) {  #backup
            $self->{'cursor'} = $self->{'linestart'};
            $$line = '';
            last;
        }
    }
}

# read >= 1 record lines terminating on matching 'pattern'.
sub scan_until_inclusive {
    my ($self, $pattern) = @_;
    #warn "scan_until_inclusive: /$pattern/\n";

    my $line  = \$self->{'line'};

    $self->{'blockstart'} = $self->{'linestart'};

    while ($self->_read_line) {
        last  if $$line =~ /$pattern/;
    }
}

# read >= 1 record lines terminating on matching 'pattern';
# skip optional $skip instances (default 0) before terminating.
sub scan_skipping_until {
    my ($self, $pattern, $skip, $key) = (@_, 0);
    #warn "scan_skipping_until: /$pattern/\n";

    my $line  = \$self->{'line'};

    $self->{'blockstart'} = $self->{'linestart'};

    while ($self->_read_line) {
        if ($$line =~ /$pattern/) {
            if ($skip < 1) {  #backup
                $self->{'cursor'} = $self->{'linestart'};
                $$line = '';
                last;
            }
            $skip--;
        }
    }
}

# read >= 1 record lines terminating on matching 'pattern'.
sub scan_while {
    my ($self, $pattern) = @_;
    #warn "scan_while: /$pattern/\n";

    my $line  = \$self->{'line'};

    $self->{'blockstart'} = $self->{'linestart'};

    while ($self->_read_line) {
        if ($$line !~ /$pattern/) {  #backup
            $self->{'cursor'} = $self->{'linestart'};
            $$line = '';
            last;
        }
    }
}

###########################################################################
# private methods
###########################################################################
#read next line of text into self.line and return 1 on success; otherwise set
#self.line to undef and return 0
sub _read_line {
    my $self = shift;
    #warn sprintf("_read_line: > linestart=%d cursor=%d extent=%d\n",
    #    $self->{'linestart'}, $self->{'cursor'}, $self->{'extent'});

    my $line = \$self->{'line'};

    $$line = undef, return 0  if $self->{'cursor'} >= $self->{'extent'};

    #read the line
    my $bytes = $self->{'text'}->getline($line, $self->{'cursor'});

    #warn "_read_line:   read $bytes\n";

    return 0  unless $bytes;

    #ignore leading indent chars
    if ($self->{'indent'} > 0) {
        $$line = substr($$line, $self->{'indent'}, $bytes - $self->{'indent'});
    }

    #advance cursor
    $self->{'linestart'} = $self->{'cursor'};
    $self->{'cursor'}   += $bytes;

    #warn sprintf("_read_line: < linestart=%d cursor=%d\n",
    #             $self->{'linestart'}, $self->{'cursor'});

    return 1;
}

###########################################################################
# debug
###########################################################################
#return the remaining unprocessed stream without altering the cursor
sub inspect_stream {
    my $self = shift;
    return ''  if $self->{'cursor'} >= $self->{'extent'};
    return $self->{'text'}->substr(
        $self->{'cursor'}, $self->{'extent'} - $self->{'cursor'} + 1
    );
}

###########################################################################
1;
