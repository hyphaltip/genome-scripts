#!/usr/bin/perl -w

# This script assumes very specific ZFF right now, that you are passing in 
# export.ann and export.dna that were processed from fathom in SNAP

use strict;
use List::Util qw(min max);
use File::Spec;
use Getopt::Long;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Location::Split;
use Bio::Location::Simple;

my $zff = shift @ARGV || 'export.ann';
my $genome = shift @ARGV || 'export.dna';

if( ! defined $zff || ! defined $genome ) {
    die "usage: zff2augustus_gbk.pl <zff> <genomefasta>\n";
}

my $dbh = Bio::DB::Fasta->new($genome);

my $out = Bio::SeqIO->new(-format => 'genbank',-fh => \*STDOUT);
open(my $fh => $zff) || die $!;

my $seq;
my @location;
while (<$fh>) {
    if (/^>(\S+)/) {
	if( $seq ) {
	    # process the previous sequence first
	    my $loc;
	    if( @location == 1 ){ 
		$loc = shift @location;
	    } else { 
		$loc = Bio::Location::Split->new();	    
		for my $locsub ( sort { $a->start <=> $b->start } @location) {
		    $loc->add_sub_Location( $locsub ); 
		}
	    }
	    my $gene = Bio::SeqFeature::Generic->new(-primary_tag => 'CDS',
		-location => $loc);
	    $seq->add_SeqFeature($gene);
	    $out->write_seq($seq);
	    @location = ();			
	}
	my $seqid = $1;
	my $seqstr = $dbh->seq($seqid);
	if( ! defined $seqstr ) {
	    die("cannot find $seqid in the input file $genome\n");
	}
	$seq= Bio::Seq->new(-seq => $seqstr, -id  => $seqid);
	$seq->add_SeqFeature(Bio::SeqFeature::Generic->new(-primary_tag => 'source',
							   -start => 1,
							   -end   => $dbh->length($seqid)));
    } else {
	my @f = split;
	if( @f != 4 && @f != 9 ) {
	    die "input does not appear to be ZFF";
	}
	my $strand = $f[1] < $f[2] ? '1' : '-1';
	if ($strand < 0) {($f[1], $f[2]) = ($f[2], $f[1])}
	my $id = pop @f;
#	warn("start,end are $f[1], $f[2]\n");
	push @location, Bio::Location::Simple->new
	    (-start => $f[1],
	     -end   => $f[2],
	     -strand => $strand);
#	warn("loc is ",$location[-1]->to_FTstring(), "\n");
    }
}

# fencepost, gotta do this again for the last one
my $loc;
if( @location == 1 ){ 
    $loc = shift @location;
} else { 
    $loc = Bio::Location::Split->new();	    
    for my $locsub ( sort { $a->start <=> $b->start } @location) {
	$loc->add_sub_Location( $locsub ); 
    }
}
my $gene = Bio::SeqFeature::Generic->new(-primary_tag => 'CDS',
					 -location => $loc);	    
$seq->add_SeqFeature($gene);
$out->write_seq($seq);

