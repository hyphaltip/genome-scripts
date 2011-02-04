#!/usr/bin/perl -w
use strict;

# create seq file
use Bio::DB::Fasta;

my $ws = 100;
my $readsize = 43;
my $db = Bio::DB::Fasta->new(shift @ARGV);
my $outname = shift || "mappable_reads.fas";
open(my $fh => ">$outname") || die $!;

for my $chrom ( $db->ids ) {
    my $len = $db->length($chrom);
    for(my $start =1; $start < $len - $readsize; $start ++ ){ 
	printf $fh ">%s.%d\n%s\n",$chrom,$start,
	$db->seq($chrom, $start => $start + $readsize-1);
    }
}
