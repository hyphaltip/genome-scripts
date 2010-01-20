#!/usr/bin/perl -w
use strict;

open(my $fh => $ARGV[0]) || die $!;
my $count = 0;
open(my $ofh => ">$ARGV[0].$count") || die $!;
my $first = <$fh>;
print $ofh $first;

while(<$fh>) {
    if(/^NSEQS/) {
	$count++;
	close($ofh);
	open($ofh => ">$ARGV[0].$count") || die $!;
    }
    print $ofh $_;    
}
