#!/usr/bin/perl
#this script should print lines that contain info for sRNA reads
#that are present in the mRNA dat file but not in the exon dat file


use strict;
use warnings;


my $mRNA_dat = $ARGV[0]; 
my $exon_dat = $ARGV[1];
my $line_count = 1;
my %exon_hash;

open (EXON, $exon_dat);
while (my $line=<EXON>){
    chomp $line;
    if ($line_count==1){
	print $line,"\n";
    }
    else{
	my @exonLine_array = split(/\t/,$line);
	$exonLine_array[8]=~s/\.exon\S//;
	$exon_hash{$exonLine_array[0]} = $exonLine_array[8];
    }
    $line_count++;
}
open (MRNA, $mRNA_dat);
while (my $line=<MRNA>){
    chomp $line;
    my @mrnaLine_array = split(/\t/,$line);
    unless(exists $exon_hash{$mrnaLine_array[0]}){
	print $line,"\n";
    }
}
