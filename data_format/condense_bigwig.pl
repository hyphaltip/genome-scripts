#!/usr/bin/perl -w
use strict;
use warnings;


=head1 NAME

condense_bigwig - summarize BEDGraph file to reduced number of windows

=head1 SYNOPSIS

condense_bigwig file.bedgraph -window 1000 > file.condense.bedgraph

=cut

use Bio::DB::BigWig 'binMean';;
use Getopt::Long;

my $window = 1000;

GetOptions(
	   'w|window:i' => \$window,
);
my ($wigfile,$genome) = @ARGV;

my $wig  = Bio::DB::BigWig->new(-bigwig=>$wigfile,
				-fasta=>$genome);

# fix this to just get windows of a fixed length?
my @bins = $wig->features(-type=>"bin:$window"); # same as "bin:1"
for my $b (@bins) {
  my $chrom = $b->seq_id;
  print join("\t",($chrom,$b->start, $b->end,$b->mean),"\n";
}
