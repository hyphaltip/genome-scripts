#!/usr/bin/perl -w
use warnings;
use strict;

=head1 NAME

assign_productnames_sprot - assign product names from Swissprot hits

=head1 DESCRIPTION

Assign product names for annotation for GenBank submission

Criteria based on Broad Institute's guidelines.
http://www.broadinstitute.org/annotation/genome/Geomyces_destructans/GeneFinding.html

"gene product names via blast hits to SwissProt database with manually
curated gene product names, using a set of strigent criteria (>=70%
protein sequence identity, >=70% coverage of both the query and the
database hit sequence, and length difference <=10%."

=head1 AUTHOR

Jason Stajich jason.stajich[at]ucr.edu

=cut

use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $url = 'http://www.uniprot.org/uniprot/%s.txt';
my $cache_dir = 'sprot_cache';

my $align_percent    = 60;
my $identity_percent = 60;
my $lendiff_percent  = 15;
my $db_file;
my $outfile;
GetOptions
    ('cache:s' => \$cache_dir,
     'url:s'   => \$url,
     'db:s'    => \$db_file, 
     'o|out:s' => \$outfile,
    );

my $sprot_tabfile = shift || die $!;
mkdir($cache_dir) unless -d $cache_dir;

open(my $fh => $sprot_tabfile) || die $!;
if( ! $outfile ) {
    if( $sprot_tabfile =~ /(\S+)\.tab$/) {
	$outfile = $1;
    } else { $outfile = $sprot_tabfile; }
    $outfile .= ".prod_names";
}
open(my $ofh => ">$outfile") || die "cannot open $outfile: $!";

my $dbh;
if( $db_file ) {
    $dbh = Bio::DB::Fasta->new($db_file);
}
my ($qname,$query_len);
my %results;
while(<$fh>) {
    if( /^\#\s+FASTA/ ) {
	$qname = $query_len = undef;
    } elsif( /^\#\s+Query:\s+(\S+).+ - (\d+) aa/) {
	$qname = $1;
	$query_len = $2;
	$results{$qname} = undef;
    } elsif ( ! $dbh && 
	      /^\#\s+Database:\s+(\S+)/ ) {
	$db_file = $1;
	$dbh = Bio::DB::Fasta->new($db_file);
    } elsif ( /^\#/ ) {
	next;
    } else {
	my ($query,$target,$pid,$similar,$mismatch,$gaps,
	    $qstart,$qend,$tstart,$tend, $evalue, $bits) = split(/\s+/,$_);
	next if $results{$query};
	next if $pid < $identity_percent;
	if( ! defined $query_len ) {
	    warn("no query len for $query ($qname)\n");
	}
	my $tlen = $dbh->length($target);
	if( 100 * (abs($tlen - $query_len) / $query_len) > $lendiff_percent ) {
	    warn("lengths are different by > 10% $query=$query_len $target=$tlen\n");
	    next;
	}
	if( 100 * ($tend - $tstart) / $tlen < $align_percent ) {
	    warn "Target $target alignment $tstart..$tend ",($tend - $tstart), " aa is less than $align_percent\n";
	    next;
	}
	if( 100 * ($qend - $qstart) / $query_len < $align_percent ) {
	    warn "Query $query alignment $qstart..$qend ",($qend - $qstart), " aa is less than $align_percent\n";
	    next;
	}
	if( $target =~ /^(sp|emb|tr)\|([^\|]+)\|(\S+)/ ) {
	    
	    my ($accession,$sprot_name) = ($2,$3);
	    if( ! $accession ) {
		die("accession not defined in $target\n");
	    }
	    my $filecache = File::Spec->catfile($cache_dir,
						$accession.".txt");
	    
	    if( ! -f $filecache ) {
		my $geturl = sprintf($url,$accession);
		`curl -L -s -o $filecache $geturl`;
	    }
	    my $seqin = Bio::SeqIO->new(-format => 'swiss',
					-file   => $filecache);
	    if( my $seq = $seqin->next_seq) {
		my $desc= $seq->description;
		my ($product_name,$ec) = ('','');
		if( $desc =~ /RecName:\s+([^;]+);/ ) {		    
		    $product_name = $1;
		    unless( $product_name =~ s/Full=// ) {
			warn("unexpected prod name '$product_name'\n");
		    }
		} 
#		$product_name =~ s/homolog(\s+\d+)?/protein/;
#		$product_name =~ s/(\s+)gene/$1protein/;
		# ad hoc stuff
		#$product_name =~ s/protein homolog//;
		if( $desc =~ /EC=([^;]+);/ ) {
		    $ec = $1;
		}

		$results{$query} = [ $target, $pid, $evalue, $product_name ];
		print $ofh join("\t", $query,$target, $pid, $evalue,$product_name,$ec), "\n";
	    }
	}
    }
}
