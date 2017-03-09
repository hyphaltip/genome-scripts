#!/usr/bin/perl
use strict;
use warnings;

if( @ARGV != 3 ) {
    die"usage: perl split_fastq_interleaved.pl INPUT.fq OUT1.fq OUT2.fq\n";
}

my ($in,$out1,$out2) = @ARGV;

if( ! -f $in ) {
    die"need an input file, you provided '$in' but it does not exist\n";
}
if( -f $out1 || -f $out2 ) {
    die "will not overwrite $out1 or $out2, please remove or provide different output file names";
}

open(my $fh => $in)|| die "cannot open $in: $!";
open(my $ofh1 => ">$out1") || die "cannot open $out1: $!";
open(my $ofh2 => ">$out2") || die "cannot open $out2: $!";

my (@ofh)  = ( $ofh1,$ofh2);
my $i =0;
while (<$fh>) {
    my $ofh = $ofh[$i];
    print $ofh $_;
    for ( 0..2 ) {
	my $line = <$fh>;
	print $ofh $line;
    }
    $i = 1 - $i; # flip flop
}
