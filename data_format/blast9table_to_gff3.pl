#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $cutoff = 1e-3;
GetOptions('c|cutoff:s' => \$cutoff);
my $i =0;
while(<>) {
    next if /^\#/;
    my @row = split;
    next if $row[10] > $cutoff;
    my ($start,$end) = sort { $a <=> $b } ($row[8],$row[9]);
    print join("\t", $row[1],'BLAST','match',
	       $start,$end,$row[2],$row[8] < $row[9] ? '+' : '-',
	       '.',
	       sprintf("ID=match%d;Name=%s;Target=%s+%d+%d",
		       $i++,$row[0],$row[0],$row[6],$row[7])),"\n";
}
    
