#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $nlen = 50;
my $in = Bio::SeqIO->new(-format => 'fasta',
						 -file   => shift @ARGV);
my @lengths;
my $total_length = 0;
my $n_total;
my $n_runs;
my $n_run_bp;
while( my $seq = $in->next_seq ) {
	my $seqstr = $seq->seq;
	$n_total += ( $seqstr =~ tr/Nn/Nn/);
	while( $seqstr =~ /([Nn]{$nlen,})/g ) {
	 $n_runs ++;
	 $n_run_bp += length($1);
	}
	push @lengths, $seq->length;
	$total_length += $seq->length;
}

print "Mean size of contigs is: ", mean_size(@lengths),"\n";
print "Median size of contigs is: ", median_size(@lengths),"\n";
print "N50 is ", find_N50(@lengths), "\n";

print "Total N bases is $n_total\n";
print "Total number of N runs >= $nlen is $n_runs for $n_run_bp\n";


sub mean_size {
	my @numbers = @_;
	my $total = 0;	
	for my $n ( @numbers ){
		$total += $n;
	}
	return $total / scalar @numbers;
}

sub median_size {
	my @numbers = @_;
	my @sorted = sort { $a <=> $b } @numbers;
	my $middle = int (scalar @sorted / 2);
	return $sorted[$middle];
}

sub find_N50 {
	my (@numbers)  = @_;
	my $total = 0;
	for my $n ( @numbers ){
		$total += $n;
	}
	my @sorted = sort { $b <=> $a } @numbers;
	my $breakpoint = 0.5 * $total;
	my $accumulate = 0;
	for my $n ( @sorted ) {
		$accumulate += $n;
		if( $accumulate >= $breakpoint) {
			return $n;
		}
	}
}
