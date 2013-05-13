#!/usr/bin/perl -w
use strict;
use warnings;


=head1 NAME

condense_bigwig - summarize BigWig file to reduced number of windows

=head1 SYNOPSIS

condense_bigwig file.bw genome.fa -window 1000 > file.condense.bedgraph

=head1 DESCRIPTION

Used to downsample for Circos drawing

=head1 AUTHOR - Jason Stajich

Jason Stajich, jason.stajich[at]ucr.edu

=cut

use Bio::DB::BigWig 'binMean';
use Getopt::Long;

my $windowsize = 1000;

GetOptions(
	   'w|windowsize:i' => \$windowsize,
);
my ($wigfile,$genome) = @ARGV;

my $wig  = Bio::DB::BigWig->new(-bigwig=>$wigfile,
				-fasta=>$genome);
my @seq_ids = $wig->seq_ids;
for my $sum ( $wig->features(-type=>'summary') ) { # one for each chromosome
    warn "chromosome ",$sum->seq_id,"\n";
    my $stats = $sum->statistical_summary($windowsize);
    warn "\tmax  = ",$stats->[0]{maxVal},"\n";
    warn "\tmean = ",binMean($stats->[0]),"\n";
}

for my $chrom ( @seq_ids ) {
    my $len = $wig->length($chrom);
    my $bincount =int( $len / $windowsize);
    warn("bincount = $bincount for $chrom:1-$len\n");

    my @bins = $wig->features(
	-seq_id => $chrom,
	-type=>"bin:$bincount");

    for my $b (@bins) {
	print join("\t",$chrom,$b->start, $b->end,$b->mean),"\n";
    }
}
