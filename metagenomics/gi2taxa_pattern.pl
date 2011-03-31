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

# the third file that is needed is GI to TaxonID
# ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
my $gi_taxid_file = 'gi_taxid_nucl_all.dmp';

my $input;
my $dbfile = '/tmp/gi2taxonomy.db';
my $force = 0;
my $cache_idx = '/tmp/taxdb';
my $debug = 0;
GetOptions('t|taxonomy:s' => \$taxonomy_folder,
	   'g|gi:s'       => \$gi_taxid_file,
	   'i|input:s'    => \$input,
	   'dbfile:s'     => \$dbfile,
	   'idx:s'        => \$cache_idx,
	   'force!'       => \$force,
	   'v|debug!'     => \$debug,
	   );
mkdir($cache_idx) unless -d $cache_idx;

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
$dbh->do(<<SQL
CREATE INDEX IF NOT EXISTS i_gi2taxonid_gi ON gi2taxonid (gi);
SQL
	 );

    $dbh->commit;
    $dbh->do("DELETE FROM gi2taxonid");
    $dbh->commit;

    my $in = File::Spec->catfile($taxonomy_folder,
				$gi_taxid_file);
    my $gifh;
    if( $in =~ /\.gz/) {
	open($gifh => "zcat $in |") || die "$in: $!\n";
    } else {
	open($gifh => "< $in") || die "$in: $!\n";
    }
    my $insert = $dbh->prepare(<<SQL
INSERT INTO gi2taxonid (gi,taxonid) VALUES (?,?)
SQL
); 
    my $count = 0;
    while(<$gifh>) {
	my ($gi,$taxid) = split;
	$insert->execute($gi,$taxid);
	$dbh->commit, warn("$count entries\n") if (++$count % $commit_interval == 0);
    }	
    $dbh->commit;    
}

warn("processing IDs\n") if $debug;
open(my $fh => $input) || die "cannot open $input: $!\n";
my $i;


my $query = $dbh->prepare(<<SQL
SELECT taxonid FROM gi2taxonid g WHERE g.gi = ?
SQL
			    );
my $missing = 0;
my (%rank,%genus,%phyla);
my (%name2id,%genus2phyla);
while( <$fh>) {
    my ($gi) = split;    
    $query->execute($gi);

    my $res = $query->fetch;
    warn("done query for gi=$gi\n") if $debug;    
    if( $res && @$res ) {
	my ($taxon) = $tdb->get_taxon(-taxonid => $res->[0]);
	if( ! defined $taxon ) {
	    warn("Cannot find taxa for taxid ", $res->[0], "\n");
	    $missing++;
	} else {
	    $rank{$taxon->rank}++;
	    my @lineage = $tree_functions->get_lineage_nodes($taxon);
	    my ($phylum) = grep { $_->rank eq 'phylum' } @lineage;
	    if( $phylum ) {
		$phyla{$phylum->scientific_name}++;
		$name2id{$phylum->scientific_name} = $phylum->id;
	    }
	    my ($g) = grep { $_->rank eq 'genus' } @lineage;
	    if( $g) {
		$genus{$g->scientific_name}++;
		$genus2phyla{$g->scientific_name} = $phylum ? $phylum->scientific_name : 'Unranked';
		$name2id{$g->scientific_name} = $g->id;
	    }
	}
    } else {	
	$missing++;
	warn("Cannot find $gi in gi2taxid\n");
    }
    $query->finish;
    last if $i++ > 10000 && $debug;
}

$dbh->disconnect;
warn("missing $missing\n");
open(my $ofh => ">rank.csv") || die $!;
print $ofh join(",",qw(RANK COUNT)),"\n";
for my $r ( sort { $rank{$b} <=> $rank{$a} } keys %rank ) {
    print $ofh join(",", $r, $rank{$r}), "\n";
}

open($ofh => ">genus.csv") || die $!;
print $ofh join(",",qw(TAXID GENUS COUNT PHYLA)),"\n";
for my $r ( sort { $genus{$b} <=> $genus{$a} } keys %genus ) {
    print $ofh join(",", $name2id{$r}, $r, $genus{$r},$genus2phyla{$r}), "\n";
}


open($ofh => ">phyla.csv") || die $!;
print $ofh join(",",qw(TAXID PHYLA COUNT)),"\n";
for my $r ( sort { $phyla{$b} <=> $phyla{$a} } keys %phyla ) {
    print $ofh join(",", $name2id{$r}, $r, $phyla{$r}), "\n";
}
