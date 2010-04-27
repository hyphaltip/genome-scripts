#!/usr/bin/perl -w
use Bio::SeqIO;
use strict;
use Getopt::Long;

# parameters
# -g or --genome to give the genome as a FASTA file
my $genome;
# -i or --input to get the name of the input file
my $input_pileup; # name of the input file in PILEUP format
my $bin_size = 1000;
my $min_read_depth = 10;
GetOptions(
	   'g|genome:s'   => \$genome,
	   'i|input:s'    => \$input_pileup,
	   'b|binsize:i'  => \$bin_size,
	   'm|mindepth:i' => \$min_read_depth,
	   );

if( ! defined $genome ) {
    die("must provide a genome file in FASTA format\n");
}

# this means you can do this
# perl heterozygote_estimate.pl -g genome inputfile
# or this
# perl heterozygote_estimate.pl -g genome -i inputfile
$input_pileup = shift @ARGV unless defined $input_pileup;
if( ! defined $input_pileup ) {
    die("must provide a pileup file\n");
}

my %chrom_lengths;
my %bins;
my $in = Bio::SeqIO->new(-format => 'fasta',
			 -file   => $genome);
# read in all the seqs one-at-a-time
while( my $seq = $in->next_seq ) {
    # store the chrom lengths, indexing by chromosome ID
    $chrom_lengths{$seq->id} = $seq->length;

    # auto-vivify (initialize the length of the array to the number of bins
    # it takes to cover each chromosome
    $bins{$seq->id}->[ int($seq->length / $bin_size)] = 0;
}

# store some data

open(my $fh => $input_pileup) || die "cannot open $input_pileup: $!";
while(<$fh>) {
    chomp; # remove the trailing \n
    my ($chrom_name, $variant_pos, $ref_allele, $variant_alleles,
	$consensus_score, $variant_score,
	$mapping_quality, $read_depth) = split(/\t/,$_);
    
#    next if( $read_depth < $min_read_depth );
    my $bin = int ( $variant_pos / $bin_size);
#    print "variant is at $variant_pos on chrom $chrom_name, goes in bin $bin\n";
    $bins{$chrom_name}->[ $bin ]++;
#    last;
}

for my $chrom ( keys %bins ) {
    for my $bin ( @{ $bins{$chrom} } ) {
	printf "%d\n", $bin || 0;  #print the number of SNPs seen
#	printf "%.4f\n", ($bin || 0) / $bin_size; # print SNP rate
    }
}
