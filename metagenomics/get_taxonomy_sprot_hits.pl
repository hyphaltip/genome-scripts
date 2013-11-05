#!/usr/bin/perl -w
use strict;


use Bio::Index::Swissprot;
use Bio::DB::SwissProt;
use Bio::Tree::Tree;
use Bio::DB::Taxonomy;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;
use File::Path qw(make_path);
use Env qw(USER);
use DB_File;
#use LWP;

#my $sprot_url = 'http://www.uniprot.org/uniprot/';
my $sprot_db;
my $cache_base = "/scratch/$USER";
my $taxa_cache = "taxdb";
my $sprot_files_cache = "sprotfiles";

my $debug = 0;
my $sp_taxa = "sprot2taxonomy.db";
my $taxonomy_folder = '/srv/projects/db/taxonomy/ncbi';

my $tree_functions = Bio::Tree::Tree->new(); # for taxonomy queries
my $query_peps;
my $outdir = "sprot_tax_groups";
GetOptions(
    'v|debug!' => \$debug,
    'q|query:s'  => \$query_peps,
    'db|sprot:s' => \$sprot_db,
    'c|cache:s'  => \$cache_base,
    't|taxonomy:s' => \$taxonomy_folder,
    'o|out:s'    => \$outdir,
    );
make_path($cache_base);
make_path($outdir);

$taxa_cache = File::Spec->catdir($cache_base,$taxa_cache);
$sprot_files_cache = File::Spec->catdir($cache_base,$sprot_files_cache);
$sp_taxa = File::Spec->catfile($cache_base,$sp_taxa);

mkdir($taxa_cache);
mkdir($sprot_files_cache);


if( ! $query_peps ) {
    die "need a query pepfile\n";
}
my $qdb = Bio::DB::Fasta->new($query_peps);
my $sdb;
my $sdb_remote = Bio::DB::SwissProt->new;
if( ! $sprot_db ) {
    warn("must provide path to Swissprot DB\n");
    $sdb = $sdb_remote;
} else {
    my $sprot_index = File::Spec->catfile($cache_base,"sprot.idx");

    if( ! -f $sprot_index ) {
	$sdb = Bio::Index::Swissprot->new(-filename => $sprot_index,
					  -write_flag => 1);
	$sdb->make_index($sprot_db);
    } else {
	$sdb = Bio::Index::Swissprot->new(-filename => $sprot_index,
					  -write_flag => 0);
    }
}

my $tdb = Bio::DB::Taxonomy->new
    (-source => 'flatfile',
     -nodesfile => File::Spec->catfile($taxonomy_folder, 'nodes.dmp'),
     -namesfile => File::Spec->catfile($taxonomy_folder, 'names.dmp'),
     -directory => $taxa_cache,
     );
my %groups;
for my $infile ( @ARGV ) {
    open(my $fh => $infile) || die "Cannot open $infile: $!";
    my $i = 0;
    my %seen;
    while(<$fh>) {
	next if /^\#/;
	my ($q,$s,$pid,@rest) = split;
	next if $seen{$q}++; # only look at the 1st hit for now
	my $qseq = $qdb->get_Seq_by_acc($q);
	my (undef,$acc,$spid) = split(/\|/,$s);
	my $bitscore = pop @rest;
	my $evalue   = pop @rest;
	
	my $seq = $sdb->get_Seq_by_id($spid);
	if( ! $seq ) {
	    # try accession
	    $seq = $sdb->get_Seq_by_acc($acc);
	    eval { 
		if( ! $seq ) {
		    # go remote
		    $seq = $sdb_remote->get_Seq_by_id($spid);
		}
	    };
	    if( $@ ) {
		warn("cannot find $s in the dbs\n");
	    }
	}
	if( $seq ) {
	    my $taxon = $seq->species;
	    my ($species_name,$genus) = ($taxon->classification);
	    my @taxonids = $tdb->get_taxonids($species_name);
	    if( ! @taxonids ) {
		@taxonids = $tdb->get_taxonids($genus);
	    } 
	    if( ! @taxonids ) {
		warn("cannot get taxonid for $species_name ($genus)\n");
	    }
	   
	    if( @taxonids > 1 ) {
		warn("multiple IDs found for $genus $species_name\n");
	    }
	    for my $t ( @taxonids ) {
		my $str = &get_rankstring($t);
		push @{$groups{$str}}, [$q,$spid,$pid,$evalue,$qseq];
	    }
	} else {
	    warn("cannot find $spid ($acc)\n");
	    next;
	}
	last if $debug && $i++ > 100;
    }
}
for my $grp ( keys %groups ) {
	my $outseq = Bio::SeqIO->new(-format => 'fasta',
				     -file   => ">$outdir/$grp.query_hits.fas");
	open(my $ofh => ">$outdir/$grp.query_hits.report") || die $!;
	for my $hit ( @{$groups{$grp}} ) {
	    my $qseq = pop @$hit;
	    $outseq->write_seq($qseq);
	    print $ofh join("\t",@$hit),"\n";
	}
	close($ofh);
}

sub get_rankstring {
    my $taxid = shift;
    
    my ($taxon) = $tdb->get_taxon(-taxonid => $taxid);

    my @lineage = $tree_functions->get_lineage_nodes($taxon);
    shift @lineage; # drop root
    my @ordered_names;
    for my $l ( @lineage ) {
	if( $l->rank =~ /^(superkingdom|kingdom|phylum)/i ) {
	    push @ordered_names, $l->scientific_name;
	}
    }
    return join("__", @ordered_names);
}
