#!/usr/bin/perl -w
use strict;
my $i = 0;
my %nums;
while(<>) {
    next unless /^\+/;
    my $line= <>;
    chomp($line);
    while( my $n = substr($line,0,1,'') ) {
	$nums{ord($n)}++;
    }
    last if $i++ > 10_000;
}
for my $n ( sort { $a <=> $b } keys %nums ) {
    printf "%d\t%d\t%s\t%s\n",$n-64,$n,chr($n),$nums{$n};
}
