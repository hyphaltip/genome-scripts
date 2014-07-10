#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use Getopt::Long;

my $debug = 0;
my $in_map = 'better.map';
my $in_genomes = 'genomes';
my $out_dir = 'circos/data';

GetOptions(
    'v|debug!'   => \$debug,
    'm|map:s'    => \$in_map,
    'g|genome:s' => \$in_genomes,
    'o|out:s'    => \$out_dir,
    );

make_path($out_dir) unless -d $out_dir;
open(my $fh => $in_genomes) || die "cannot open $in_genomes: $!";

my @genomes;
while(<$fh>) {
    @genomes = split;
    last;
}
close($fh);
warn join("\n",@genomes),"\n" if $debug;

my @map;
open($fh => $in_map) || die "cannot open $in_map: $!";
while(<$fh>) {
    my ($block_id, @sets) = split;
    for my $genome ( @genomes ) {
	my ($contig,$start,$end,$strand) = splice(@sets, 0,4);
	warn("contig is $contig\n") if $debug;
	next if $contig eq 'NA';
	$contig = $genome. "_" . $contig;
	$map[$block_id]->{$genome} = [ $contig, $start, $end, $strand ];
    }
}

# output all pairwise block information
for( my $i = 0; $i < @genomes; $i++ ) {
    my $pairA = $genomes[$i];
    for( my $j = $i+1; $j < @genomes; $j++ ) {
	next if $i == $j;
	my $pairB = $genomes[$j];
	open(my $ofh => ">$out_dir/$pairA-$pairB.all.txt") || die $!;
	warn "$pairA-$pairB\n" if $debug;
	my $n = 1;
	for my $block ( @map ) {
	    next if ! defined $block;
	    if( exists $block->{$pairA} &&
		exists $block->{$pairB} ) {
		print $ofh join(" ", 
				(map { $block->{$pairA}->[$_] } (0,1,2) ), # chr start stop
				(map { $block->{$pairB}->[$_] } (0,1,2) ),
				sprintf("id=%d",$n)
		    ), "\n";				
	    }
	    $n++;
	}
    }
}
    
