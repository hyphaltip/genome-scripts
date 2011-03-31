#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::Taxonomy;
use File::Spec;
use Bio::DB::GenBank;
use Bio::Tree::Tree;
use DBI qw(:sql_types);

my $remote = Bio::DB::GenBank->new;
my $tree_functions = Bio::Tree::Tree->new(); # for taxonomy queries
my $commit_interval = 500_000;

# Author Jason Stajich - jason[at]bioperl.org
# This will take a FASTA file (or just FASTA headers) with 
# embedded accession numbers (as first part of FASTA header, before the '_'
# and converts to a taxonomy file suitable for QIIME, e.g.:
# AY800210 Archae; Euryarchaeota; Halobacteriales


# this implementation currently uses SQLite

my $taxonomy_folder = '/project/db/taxonomy/ncbi';
# This folder needs to contain the uncompressed contents of 
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# The livelist (accession numbers to GI numbers) needs to be stored in this
# folder as well.
# I is ftp://ftp.ncbi.nih.gov/genbank/livelists/GbAccList.1003.2010.gz'
# so you can do wget or lftp for this file
my $livelist_file = 'GbAccList.1003.2010.gz';

# the third file that is needed is GI to TaxonID
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
my $gi_taxid_file = 'gi_taxid_nucl.dmp';

my $input;
my $dbfile = '/tmp/blastids2taxonomy.db';
my $force = 0;
my $cache_idx = '/tmp';
my $debug = 0;
GetOptions('t|taxonomy:s' => \$taxonomy_folder,
	   'l|livelist:s' => \$livelist_file,
	   'g|gi:s'       => \$gi_taxid_file,
	   'i|input:s'    => \$input,
	   'dbfile:s'     => \$dbfile,
	   'idx:s'        => \$cache_idx,
	   'force!'       => \$force,
	   'v|debug!'     => \$debug,
	   );

$input ||= shift @ARGV;
die "need an input file with -i \n" unless defined $input;

$force = 1 unless( -f $dbfile ); # whether or not to remake the lookup tables in SQL

my $tdb = Bio::DB::Taxonomy->new
    (-source => 'flatfile',
     -nodesfile => File::Spec->catfile($taxonomy_folder, 'nodes.dmp'),
     -namesfile => File::Spec->catfile($taxonomy_folder, 'names.dmp'),
     -directory => $cache_idx,
     );

warn("initializing DBH\n") if $debug;
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","", 
		       {AutoCommit => 0});
warn("done initializing DBH\n") if $debug;


if( $force ) {
    warn("rebuilding DB\n");
    $dbh->do(<<SQL
CREATE TABLE IF NOT EXISTS gi2taxonid ( 
 gi2taxonid_id  INTEGER PRIMARY KEY ASC, 
 gi      INTEGER NOT NULL,
 taxonid INTEGER NOT NULL
);
	 SQL
	     );
    $dbh->do(<<SQL
CREATE UNIQUE INDEX IF NOT EXISTS ui_gi2taxonid ON gi2taxonid (gi,taxonid);
SQL
	     );
#$dbh->do(<<SQL
#CREATE INDEX IF NOT EXISTS i_gi2taxonid_gi ON gi2taxonid (gi);
#SQL
#	 );

    $dbh->do(<<SQL
CREATE TABLE IF NOT EXISTS acc2gi (
 acc2gi_id INTEGER PRIMARY KEY ASC,
 accession varchar(24) NOT NULL,
 gi        INTEGER NOT NULL
				   );
SQL
	     );
    $dbh->do(<<SQL
CREATE UNIQUE INDEX IF NOT EXISTS ui_acc2gi ON acc2gi (accession);
SQL
	     );
    $dbh->do(<<SQL
CREATE INDEX IF NOT EXISTS i_gi ON acc2gi (gi);
SQL
	     );
    $dbh->commit;
    $dbh->do("DELETE FROM gi2taxonid");
    $dbh->do("DELETE FROM acc2gi");
    $dbh->commit;
    my $in = File::Spec->catfile($taxonomy_folder,
				 $livelist_file);
    my $accfh;
    if( $in =~ /\.gz/) {
	open($accfh => "zcat $in |") || die "$in: $!\n";
    } else {
	open($accfh => "< $in") || die "$in: $!\n";
    }
    $dbh->commit;
    my $insert = $dbh->prepare(<<SQL
INSERT INTO acc2gi (accession,gi) VALUES (?,?)
SQL
); 
    my ($acc,$gi,$taxid);
    my $count = 0;
    while(<$accfh>) {
	chomp;
	($acc,undef,$gi) = split(/,/,$_);
	$insert->execute($acc,$gi);
	$dbh->commit, warn("$count entries\n") if (++$count % $commit_interval == 0);
    }	
    $dbh->commit;

    $in = File::Spec->catfile($taxonomy_folder,
				$gi_taxid_file);
    my $gifh;
    if( $in =~ /\.gz/) {
	open($gifh => "zcat $in |") || die "$in: $!\n";
    } else {
	open($gifh => "< $in") || die "$in: $!\n";
    }
    $insert = $dbh->prepare(<<SQL
INSERT INTO gi2taxonid (gi,taxonid) VALUES (?,?)
SQL
); 
    $count = 0;
    while(<$gifh>) {
	($gi,$taxid) = split;
	$insert->execute($gi,$taxid);
	$dbh->commit, warn("$count entries\n") if (++$count % $commit_interval == 0);
    }	
    $dbh->commit;    
}

warn("processing IDs\n") if $debug;
open(my $fh => $input) || die "cannot open $input: $!\n";
my $i;
my $query = $dbh->prepare(<<SQL
SELECT taxonid FROM gi2taxonid g, acc2gi a
WHERE a.accession = ? AND g.gi = a.gi
SQL
			  );
my $querygi = $dbh->prepare(<<SQL
SELECT taxonid FROM gi2taxonid g WHERE g.gi = ?
SQL
			    );

while( <$fh>) {
    next unless /^>(\S+)/;
    my $desc = $1;
    next unless ( $desc =~ /^(XM_\d+|[^_]+)_/);
    my $accession = $1;
    if( ! $accession ) {
	warn("no accession in $desc\n");
	next;
    }
    warn("querying for $accession\n") if $debug;
    $query->execute($accession);

    my $res = $query->fetch;
    warn("done query for $accession\n") if $debug;    
    if( $res && @$res ) {
	&print_taxonomy($desc,$res->[0]);
    } else {	
	my $seq = $remote->get_Seq_by_acc($accession);
	if( ! $seq ) {
	    warn("no results for $accession\n");
	} else {
	    $querygi->execute($seq->primary_id); # GI is stored as primary_id in bioperl objects from GenBank parsing
	    $res = $querygi->fetch;
	    if( $res && @$res ) {	       
		&print_taxonomy($desc,$res->[0]);
	    } else {
		my @lineage = reverse $seq->species->classification;
		print join("\t", $desc, join("; ", @lineage)), "\n";
	    }
	    $querygi->finish;
	}
    }
    $query->finish;
    last if $i++ > 10 && $debug;
}

$dbh->disconnect;
	

sub print_taxonomy {
    my ($accession,$taxid) = @_;
    my ($taxon) = $tdb->get_taxon(-taxonid => $taxid);
    
    my @lineage = $tree_functions->get_lineage_nodes($taxon);
    shift @lineage; # drop root
    print join("\t", $accession, join("; ",map { $_->scientific_name } @lineage)), "\n";
}
