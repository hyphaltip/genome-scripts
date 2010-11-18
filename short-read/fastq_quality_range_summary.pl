#!/usr/bin/perl -w
use strict;
use Getopt::Long;
# print the counts of quality scores to quickly summarize quality of a FASTQ file
# default is to only process 1st 10,000 sequences
# but if --nodebug is provided as option will process all seqs in a file

my $debug = 1;
GetOptions(
	'd|debug!' => \$debug);

my $i = 0;
my %nums;
while(<>) {
    next unless /^\+/;
    my $line= <>;
    chomp($line);
    while( my $n = substr($line,0,1,'') ) {
	$nums{ord($n)}++;
    }
    last if $debug && $i++ > 10_000;
}
for my $n ( sort { $a <=> $b } keys %nums ) {
    printf "%d\t%d\t%s\t%s\n",$n-64,$n,chr($n),$nums{$n};
}
