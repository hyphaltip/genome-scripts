#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bit::Vector;
# CDS.bed 0.159%, most-cons.bed 0.233%, both 0.102%, cover 64.25%, enrich 275.73x

my %seqs;
my $ext ='exon';
my $debug = 0;
my $mapfile; # mercator mapfile
my $persite = 0;
GetOptions(
	   'm|mapfile:s' => \$mapfile,
	   'v|verbose!'  => \$debug,
	   'e|s|ext:s'   => \$ext,
	   'p|persite!'  => \$persite,
	   );
die"Need a Mapfile\n" unless defined $mapfile;

my ($bed1,$bed2)= @ARGV;
open(my $mapfh => $mapfile) || die "$mapfile: $!";
my $ref_genome_index = 0;
my @aln_lookups;
while(<$mapfh>) {
    my ($aln_id, @line) = split;
    next if $line[0] eq 'NA';
    $aln_lookups[$aln_id] = [map { $line[$ref_genome_index + $_] } 0..2 ];
}

open(my $fh => $bed1) || die $!;
my $i = 0;
my @aln_bits;
while(<$fh>) {
    next if /^\#/;
    
    next unless( !defined $ext || /\Q$ext\E/);
    my ($chrom,$start,$end,$group) = split;
    if( $chrom !~ /^\d+$/ ) {
	die("expected alignment id for chromosome\n");
    }    
    unless( exists $aln_lookups[$chrom] ){
	die("cannot find alnid '$chrom', perhaps you have the wrong Map file\n");
    }
    my $block_length = $aln_lookups[$chrom]->[2] - $aln_lookups[$chrom]->[1];

    if( ! exists $aln_bits[$chrom] ) {
	warn("block length is $block_length for '$chrom'\n") if $debug;
	$aln_bits[$chrom] = Bit::Vector->new($block_length+1);
    }
    my $offset = $aln_lookups[$chrom]->[1];
    warn(sprintf "setting %d..%d to 1 (offset is $offset alnid is '$chrom', feat is $group)\n",
	 $start,$end) if $debug;

    if( $start> $block_length) {
	next; # a feature that overlaps boundary should have been truncated
    } elsif( $end > $block_length ) {
	$end = $block_length;
    }
    $aln_bits[$chrom]->Interval_Fill($start,$end);
    $i++;
}
close($fh);
open($fh => $bed2) || die $!;
my @test_vecs;
while(<$fh>) {
    my ($chrom,$start,$end,$group) = split;
    next unless( exists $aln_bits[$chrom] ); # skip contigs with no annot features

    unless( exists $test_vecs[$chrom] ) {
	$test_vecs[$chrom] = $aln_bits[$chrom]->Shadow;
    }
    warn("Test vec: $chrom:$start..$end $group\n") if $debug;
    $test_vecs[$chrom]->Interval_Fill($start, $end);
}
my $chrom = 0;
my ($feat_sum, $test_sum,$total_len,$union_sum);
for my $vec ( @test_vecs ) {
    if( defined $vec) {
	my $vec3 = $vec->Shadow;
	$vec3->And($vec,$aln_bits[$chrom]);	
	my $feat_ct = popcount($aln_bits[$chrom]);
	$feat_sum += $feat_ct;

	my $test_ct = popcount($vec);
	$test_sum += $test_ct;

	my $len = $vec->Size;

	$total_len += $len;
	
	my $union_ct = popcount($vec3);
	$union_sum += $union_ct;

        # CDS.bed 0.159%, most-cons.bed 0.233%, both 0.102%, cover 64.25%, enrich 275.73x
	if( $persite ) {
	    printf "# AlnBlock %d, %s %5.3f%%, %s %5.3f%%, both %5.3f%%, cover %4.2f%%, enrich %4.2fx\n",
	    $chrom,
	    $bed1, 100 * ($feat_ct / $len),
	    $bed2, 100 * ($test_ct / $len),
	    100 * $union_ct / $len,
	    100 * ($union_ct / $feat_ct),
	    ($union_ct / $test_ct) / ( $feat_ct / $aln_bits[$chrom]->Size);	
	}
    }
    $chrom++;
}

# we want to print something like this:	
# CDS.bed 0.159%, most-cons.bed 0.233%, both 0.102%, cover 64.25%, enrich 275.73x
printf  "# %s %.3f%%, %s %.3f%%, both %.3f%%, cover %.2f%%, enrich %4.2fx\n",
    $bed1, 100 * ($feat_sum / $total_len),
    $bed2, 100 * ($test_sum / $total_len),
    100 * $union_sum / $total_len,
    100 * ($union_sum / $feat_sum),
    ($union_sum / $test_sum) / ( $feat_sum / $total_len);	

# from featureBits.c - this is what we want to implement
#       fprintf(stderr,"%s %5.3f%%, %s %5.3f%%, both %5.3f%%, cover %4.2f%%, enrich %4.2fx\n",
#                tables[0], 
#                100.0 * totalFirstBits/totalBases,
#                tables[1],
#                100.0 * totalSecondBits/totalBases,
#                100.0 * totalBits/totalBases,
#                100.0 * totalBits / totalFirstBits,
#                (totalBits/totalSecondBits) / (totalFirstBits/totalBases) );


# this seems fast enough
sub popcount { # hamming weight
    my $x = shift;
    my $count = 0; 
    my $size = $x->Size;
    for( my $i = 0; $i < $size; $i++ ) {
	$count++ if $x->bit_test($i);;
    }
    return $count;
}
