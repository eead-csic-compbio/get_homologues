package PlantCompUtils;
require Exporter;

# Copyright [2019-2020] EMBL-European Bioinformatics Institute
# https://github.com/Ensembl/plant_tools/tree/master/phylogenomics 

@ISA       = qw(Exporter);
@EXPORT_OK = qw(
  list_ensembl_mysql_dbs get_canonical_transcript_ids 
  download_FASTA_file parse_isoform_FASTA_file
  download_compara_TSV_file download_GTF_file get_gene_coords_GTF_file
  perform_rest_action transverse_tree_json
  minimize_MSA write_boxplot_file factorial fisher_yates_shuffle
  @DIVISIONS $REQUEST_COUNT $COMPARADIR $FASTADIR $GTFDIR $FTPURL
);

use strict;
use warnings;
use Net::FTP;
use Time::HiRes;
use HTTP::Tiny;
use DBI;

# Fungi Protists Metazoa have collections and one all-vs-all TSV file
# This code won't work there
our @DIVISIONS  = qw( Plants );
our $FTPURL     = 'ftp.ensemblgenomes.org';
our $COMPARADIR = '/pub/xxx/current/tsv/ensembl-compara/homologies';
our $FASTADIR   = '/pub/current/xxx/fasta';
our $GTFDIR     = '/pub/current/xxx/gtf';

my  $MYSQLURL   = 'mysql-eg-publicsql.ebi.ac.uk';
my  $MYSQLPORT  = 4157; 
my  $MYSQLUSER  = 'anonymous';

our $REQUEST_COUNT = 0;

my $GZIPEXE = 'gzip';  #default system gzip

# Connects to Ensembl public mysql server,
# returns reference to list of available databases
sub list_ensembl_mysql_dbs {

  my $dbh = DBI->connect("dbi:mysql:host=$MYSQLURL;port=$MYSQLPORT",$MYSQLUSER)
    || die "# ERROR: cannot connect to MySQL server $MYSQLURL port $MYSQLPORT\n";

  my $res = $dbh->selectall_arrayref('SHOW DATABASES');
  my @dbnames = map { $_->[0] } @$res;

  $dbh->disconnect();

  return \@dbnames;
}

## Connects to Ensembl public mysql server, 
## queries species core schema,
## returns reference to hash of canonical transcript stable_ids
sub get_canonical_transcript_ids {

  my ($species_core_db, $verbose) = @_;
  
  my (@tmpdata,%canon_isofs);

  my $dbh = DBI->connect("DBI:mysql:$species_core_db:host=$MYSQLURL;port=$MYSQLPORT",$MYSQLUSER)
    || die "# ERROR: cannot connect to MySQL server $MYSQLURL port $MYSQLPORT\n";

  my $sth = $dbh->prepare( 
    "SELECT t.stable_id FROM gene g, transcript t ".
    "WHERE g.canonical_transcript_id = t.transcript_id"
  );

  $sth->execute() ||
  	die "# ERROR: cannot query $MYSQLURL port $MYSQLPORT\n";

  while(@tmpdata = $sth->fetchrow_array()) {
    $canon_isofs{$tmpdata[0]} = 1;
    if($verbose){ print "$tmpdata[0]\n" }
  }

  $sth->finish();
  $dbh->disconnect();

  return \%canon_isofs;
}

# parses a FASTA file, either pep or cdna, downloaded with download_FASTA_file
# returns:
# i) a isoform=>sequence hash with the (optionally) selected or
#    (default) longest peptide/transcript per gene
# ii) a isoform=>description hash with description & chr location
#    extracted from the FASTA header
sub parse_isoform_FASTA_file {

    my ( $FASTA_filename, $ref_isoforms2keep ) = @_;

    my ( $stable_id, $gene_stable_id, $max, $len );
    my ( $iso_selected, $len_selected );
    my ( %sequence, %sequence4gene, %header, %header4gene );

    open( FASTA, "$GZIPEXE -dc $FASTA_filename |" )
      || die "# ERROR(parse_isoform_FASTA_file): cannot open $FASTA_filename\n";
    while ( my $line = <FASTA> ) {

        #>g00297.t1 pep supercontig:Ahal2.2:FJVB01000001.1:1390275:1393444:1 gene:g00297 ...
        if ( $line =~ m/^>(\S+).*?gene:(\S+)/ ) {
            $stable_id      = $1;    # might pep or cdna id
            $gene_stable_id = $2;

            if ( $line =~ m/description:(.*)?\[/ ) {
                $header{$stable_id} = $1;
            }
            else {
                $header{$stable_id} = "no_description ";
            }
            $header{$stable_id} .= ( split( /\s+/, $line ) )[2];

        }
        elsif ( $line =~ m/^(\S+)/ ) {
            $sequence{$gene_stable_id}{$stable_id} .= $1;
        }
    }
    close(FASTA);

    foreach $gene_stable_id ( keys(%sequence) ) {

        # work out which isoform should be kept for this gene
        ( $max, $iso_selected, $len_selected ) = ( 0, '', '' );
        foreach $stable_id ( keys( %{ $sequence{$gene_stable_id} } ) ) {

            # find longest isoform (default), note that key order is random
            $len = length( $sequence{$gene_stable_id}{$stable_id} );
            if ( $len > $max ) {
                $max          = $len;
                $len_selected = $stable_id;
            }

            if ( $ref_isoforms2keep->{$stable_id} ) {
                $iso_selected = $stable_id;
            }
        }

        if ($iso_selected) {
            $sequence4gene{$iso_selected} =
              $sequence{$gene_stable_id}{$iso_selected};
            $header4gene{$iso_selected} = $header{$iso_selected};
        }
        elsif ($len_selected) {
            $sequence4gene{$len_selected} =
              $sequence{$gene_stable_id}{$len_selected};
            $header4gene{$len_selected} = $header{$len_selected};
        }
        else {
            print "# ERROR(parse_isoform_FASTA_file): cannot select ".
                "an isoform for gene $gene_stable_id\n";
        }
    }

    return ( \%sequence4gene, \%header4gene );
}

# reads a compressed GTF file and returns a reference to an array with
# chr-order geneids and their coordinates
sub get_gene_coords_GTF_file {

    my ($GTF_filename) = @_;

    my ( $chr, $start, $end, $strand, $geneid );
    my (@chr_sorted_gene_ids);

    open( GTF, "$GZIPEXE -dc $GTF_filename |" )
      || die "# ERROR(get_gene_coords_GTF_file): cannot open $GTF_filename\n";
    while ( my $line = <GTF> ) {

        #1 araport11 gene 3631 5899 . + . gene_id "AT1G01010";...
        if ( $line =~
               m/^([^#])\t[^\t]+\tgene\t(\d+)\t(\d+)\t[^\t]\t(\S+)\t[^\t]\tgene_id "([^";]+)/
          )
        {
            ( $chr, $start, $end, $strand, $geneid ) = ( $1, $2, $3, $4, $5 );
            push( @chr_sorted_gene_ids,
                [ $geneid, $chr, $start, $end, $strand ] );
        }
    }
    close(GTF);

    return \@chr_sorted_gene_ids;
}

# download compressed GTF file from FTP site, renames it
# # and saves it in $targetdir; uses FTP globals defined above
sub download_GTF_file {

    my ( $dir, $ref_genome, $targetdir ) = @_;
    my ( $gtf_file, $stored_gtf_file ) = ( '', '' );

    if ( my $ftp =
        Net::FTP->new( $FTPURL, Passive => 1, Debug => 0, Timeout => 60 ) )
    {
        $ftp->login( "anonymous", '-anonymous@' )
          || die "# ERROR(download_GTF_file): cannot login " . $ftp->message();
        $ftp->cwd($dir)
          || die
          "# ERROR(download_GTF_file): cannot change working directory to $dir "
          . $ftp->message();
        $ftp->cwd($ref_genome)
          || die "# ERROR(download_GTF_file): cannot find $ref_genome in $dir "
          . $ftp->message();

        # find out which file is to be downloaded and
        # work out its final name
        foreach my $file ( $ftp->ls() ) {
            if ( $file =~ m/gtf.gz/ && $file !~ m/.chr.gtf/ ) {
                $gtf_file        = $file;
                $stored_gtf_file = "$targetdir/$gtf_file";
                last;
            }
        }

        unless ( -s $stored_gtf_file ) {
            $ftp->binary();
            my $downsize = $ftp->size($gtf_file);
            $ftp->hash( \*STDOUT, $downsize / 20 ) if ($downsize);
            printf( "# downloading %s (%1.1fMb) ...\n",
                $gtf_file, $downsize / ( 1024 * 1024 ) );
            print "# [        50%       ]\n# ";
            if ( !$ftp->get($gtf_file) ) {
                die
                  "# ERROR(download_GTF_file): failed downloading $gtf_file\n";
            }

            # rename file to final name
            rename( $gtf_file, $stored_gtf_file );

            print "# using $gtf_file\n\n";
        }
        else {
            print "# re-using $stored_gtf_file\n\n";
        }
    }
    else {
        die "# ERROR(download_GTF_file): cannot connect to $FTPURL , please try later\n";
    }

    return $stored_gtf_file;
}

# download compressed TSV file from FTP site, renames it
# and saves it in $targetdir; uses FTP globals defined above
# NOTE: if species file is not found it tries the bulky all-vs-all file
sub download_compara_TSV_file {

    my ( $dir, $ref_genome, $targetdir ) = @_;
    my ( $compara_file, $stored_compara_file ) = ( '', '' );

    if ( my $ftp =
        Net::FTP->new( $FTPURL, Passive => 1, Debug => 0, Timeout => 60 ) )
    {
        $ftp->login( "anonymous", '-anonymous@' )
          || die "# ERROR(download_compara_TSV_file): cannot login "
          . $ftp->message();
        $ftp->cwd($dir)
          || die "# ERROR(download_compara_TSV_file): cannot change working directory to $dir "
          . $ftp->message();

        # find out which file is to be downloaded
        if ( $ftp->cwd($ref_genome) ) {
            foreach my $file ( $ftp->ls() ) {
                if ( $file =~ m/protein_default.homologies.tsv.gz/ ) {
                    $compara_file        = $file;
                    $stored_compara_file = "$targetdir/$compara_file";
                    $stored_compara_file =~ s/tsv.gz/$ref_genome.tsv.gz/;
                    last;
                }
            }
        }
        else {    # try all-vs-all file instead (Fungi, Protists, Metazoa)

            print "# WARNING(download_compara_TSV_file): cannot find ".
			    "$ref_genome in $dir, try all-vs-all\n";

            foreach my $file ( $ftp->ls() ) {
                if ( $file =~ m/protein_default.homologies.tsv.gz/ ) {
                    $compara_file        = $file;
                    $stored_compara_file = "$targetdir/$compara_file";
                    foreach my $div (@DIVISIONS) {
                        if ( $dir =~ m/($div)/i ) {
                            $div = $1;
                            $stored_compara_file =~ s/tsv.gz/$div.tsv.gz/;
                            last;
                        }
                    }
                }
            }
        }

        if ( !$compara_file ) {
            die "# ERROR(download_compara_TSV_file): cannot find file to download\n";
        }

        # download that TSV file
        unless ( -s $stored_compara_file ) {
            $ftp->binary();
            my $downsize = $ftp->size($compara_file);
            $ftp->hash( \*STDOUT, $downsize / 20 ) if ($downsize);
            printf( "# downloading %s (%1.1fMb) ...\n",
                $stored_compara_file, $downsize / ( 1024 * 1024 ) );
            print "# [        50%       ]\n# ";
            if ( !$ftp->get($compara_file) ) {
                die "# ERROR(download_compara_TSV_file): failed downloading $compara_file\n";
            }

            # rename file to final name
            rename( $compara_file, $stored_compara_file );

            print "# using $stored_compara_file\n\n";
        }
        else {
            print "# re-using $stored_compara_file\n\n";
        }
    }
    else {
        die "# ERROR(download_compara_TSV_file): cannot connect to ".
            "$FTPURL , please try later\n";
    }

    return $stored_compara_file;
}

# download compressed FASTA file from FTP site, and saves it in $targetdir;
# uses FTP globals defined above
sub download_FASTA_file {

    my ( $dir, $genome_folder, $targetdir ) = @_;
    my ( $fasta_file, $stored_fasta_file ) = ( '', '' );

    if ( my $ftp =
        Net::FTP->new( $FTPURL, Passive => 1, Debug => 0, Timeout => 60 ) )
    {
        $ftp->login( "anonymous", '-anonymous@' )
          || die "# ERROR(download_FASTA_file): cannot login "
          . $ftp->message();
        $ftp->cwd($dir)
          || die "# ERROR(download_FASTA_file): cannot change working directory to $dir "
          . $ftp->message();
        $ftp->cwd($genome_folder)
          || die "# ERROR(download_FASTA_file): cannot find $genome_folder in $dir "
          . $ftp->message();

        # find out which file is to be downloaded and
        # work out its final name
        foreach my $file ( $ftp->ls() ) {
            if ( $file =~ m/all.fa.gz/ ) {
                $fasta_file        = $file;
                $stored_fasta_file = "$targetdir/$fasta_file";
                last;
            }
        }

        # download that FASTA file
        unless ( -s $stored_fasta_file ) {
            $ftp->binary();
            my $downsize = $ftp->size($fasta_file);
            $ftp->hash( \*STDOUT, $downsize / 20 ) if ($downsize);
            printf( "# downloading %s (%1.1fMb) ...\n",
                $fasta_file, $downsize / ( 1024 * 1024 ) );
            print "# [        50%       ]\n# ";
            if ( !$ftp->get($fasta_file) ) {
                die "# ERROR(download_FASTA_file): failed downloading $fasta_file\n";
            }

            # rename file to final name
            rename( $fasta_file, $stored_fasta_file );

            print "# using $fasta_file\n\n";
        }
        else {
            print "# re-using $stored_fasta_file\n\n";
        }
    }
    else {
        die "# ERROR(download_FASTA_file): cannot connect to ".
             "$FTPURL , please try later\n";
    }

    return $stored_fasta_file;
}

# uses global $REQUEST_COUNT
# takes a HTTP::Tiny object as first param
# based on examples at 
# https://github.com/Ensembl/ensembl-rest/wiki/Example-Perl-Client
sub perform_rest_action {
    my ( $http, $url, $headers ) = @_;

    $headers ||= {};
    $headers->{'Content-Type'} = 'application/json'
      unless exists $headers->{'Content-Type'};

    my ( $diff, $last_request_time, $current_time ) = ( 0, 0 );

    if ( $REQUEST_COUNT == 15 ) {    # check every 15
        $current_time = Time::HiRes::time();
        $diff         = $current_time - $last_request_time;

        # if less than a second then sleep for the remainder of the second
        if ( $diff < 1 ) {
            Time::HiRes::sleep( 1 - $diff );
        }

        # reset
        $last_request_time = Time::HiRes::time();
        $REQUEST_COUNT     = 0;
    }

    my $response = $http->get( $url, { headers => $headers } );
    my $status = $response->{status};

    if ( !$response->{success} ) {

     # check for rate limit exceeded & Retry-After (lowercase due to our client)
        if ( ( $status == 429 || $status == 599 )
            && exists $response->{headers}->{'retry-after'} )
        {
            my $retry = $response->{headers}->{'retry-after'};
            Time::HiRes::sleep($retry);

            # after sleeping see that we re-request
            return perform_rest_action( $url, $headers );
        }
        else {
            my ( $status, $reason ) =
              ( $response->{status}, $response->{reason} );
            die "# ERROR(perform_rest_action): failed REST request $url\n".
                 "# Status code: ${status}\n# Reason: ${reason}\n# Please re-run";
        }
    }

    $REQUEST_COUNT++;

    if   ( length( $response->{content} ) ) { return $response->{content} }
    else                                    { return '' }
}

# write TSV file appropriate for R boxplot function
sub write_boxplot_file {

    my ( $outfile, $n_genomes, $n_samples, $ref_data ) = @_;

    my ( $s, $sp );

    open( BOXDATA, ">", $outfile )
      || die "# ERROR(write_boxplot_file): cannot create $outfile\n";
    for ( $sp = 0 ; $sp < $n_genomes ; $sp++ ) {
        printf( BOXDATA "g%d\t", $sp + 1 );    #g = genome
    }
    print BOXDATA "\n";

    for ( $s = 0 ; $s < $n_samples ; $s++ ) {
        for ( $sp = 0 ; $sp < $n_genomes ; $sp++ ) {
            print BOXDATA "$ref_data->[$s][$sp]\t";
        }
        print BOXDATA "\n";
    }
    close(BOXDATA);

    return $outfile;
}

sub factorial {
    my $max = int( $_[0] );
    my $f   = 1;
    for ( 2 .. $max ) { $f *= $_ }
    return $f;
}

# generates a random permutation of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my ( $i, $j, $array_string );

    for ( $i = @$array ; --$i ; ) {
        $j = int( rand( $i + 1 ) );
        next if $i == $j;
        @$array[ $i, $j ] = @$array[ $j, $i ];
    }

    return join( '', @$array );
}

# takes 
# i) decoded JSON->{'tree'}, as returned from 
#          http://rest.ensembl.org/genetree/member/id/ and
# ii) reference to hash which gets populated with { accession, sequence } pairs
sub transverse_tree_json {

    my ( $tree, $ref_hash ) = @_;

    my ( $child, $id );

    if ( $tree->{'sequence'} ) {
        foreach $id ( @{ $tree->{'sequence'}{'id'} } ) {
            $ref_hash->{ $id->{'accession'} } =
              $tree->{'sequence'}{'mol_seq'}{'seq'};
        }
    }

    foreach $child ( @{ $tree->{'children'} } ) {
        transverse_tree_json( $child, $ref_hash );
    }
}

# takes a two-key hash with aligned sequences
# and returns a hash ref with minimized aligned seqs,
# after removing gap-only columns
sub minimize_MSA {

    my ( $ref_hash ) = @_;

    my $n_of_seqs = scalar( keys( %$ref_hash ));
    my ( %gap, %minMSA, $key1, $key2, $pos );

    # save position of observed gaps
    foreach $key1 (keys( %$ref_hash )) {
        foreach $key2 (keys( %{$ref_hash->{$key1}} )) {
            while($ref_hash->{$key1}{$key2} =~ m/-/g) {
                $gap{ pos($ref_hash->{$key1}{$key2})-1 }++;
            } 
        }
    }

    # copy sequences to output array but skip gap-only columns
    foreach $key1 (keys( %$ref_hash )) {
        foreach $key2 (keys( %{$ref_hash->{$key1}} )) {
            foreach $pos (0 .. length($ref_hash->{$key1}{$key2})-1) {
                next if($gap{ $pos } && $gap{ $pos } == $n_of_seqs);
                $minMSA{$key1}{$key2} .= substr($ref_hash->{$key1}{$key2},$pos,1); 
            }
        }
    }

    return \%minMSA;
}

1;
