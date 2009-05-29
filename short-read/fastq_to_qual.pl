#!/usr/bin/perl -w
use strict;
my $input = shift;
open(my $fh => $input) || die $!;
my $base = $input;
$base =~ s/\.(fq|fastq)$//;
open(my $fastafh => ">$base.fasta") || die $!;
open(my $qualfh => ">$base.fasta.qual") || die $!;
my $ready = 0;
while(<$fh>) {
    if($ready == 0 ) {
	if( /^@(\S+)/) {
	    print $fastafh ">$1\n";
	    print $qualfh  ">$1\n";
	} else {
	    die("off track $_");
	}
    } elsif( $ready == 1 ) {
	print $fastafh $_;
    } elsif( $ready == 3 ) { 
	chomp;
	my $qual = $_;
	my @row;
	while(length($qual)>0 ) {
	    my $sq = substr($qual,0,1,'');
#	    my $Q = 10 * log(1 + 10 ** (ord($sq) - 64) / 10.0) / log(10);
	    my $Q = ord($sq) - 33;
	    push @row, $Q;
	}
	print $qualfh join(" ",@row),"\n";
    } 
    $ready++;
    $ready = 0 if $ready == 4;
}
