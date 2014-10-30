#!/usr/bin/perl

# Based on original script by Peter Arensburger.
# Takes a SAM file and report location of U1A10 overlaps

# modifications by Jason Stajich
# uses flag and bitwise operator to determine +/- strand
# will read STDIN now so can pass in data on the fly
use strict;

my %topu1; #holds the location of U1 on top strand as key, name of sequence as value
my %botu1; #holds the location of U1 on bottom strand as key, name of sequence as value

#read the file, store the locations

while (my $line = <>) {
  next if $line =~ /^\@/;
  my ($name,$ori,$contig,$loc,@row) = split(/\t/,$line);
  my $seq = $row[5];
  if ($ori & 16) {
    if (substr($seq, -1, 1) eq "A") {
      $botu1{$contig}->{$loc + (length $seq) - 1} .= "$name,";
    }
  } else {
    if (substr($seq, 0, 1) eq "T") {
      $topu1{$contig}->{$loc} .= "$name,";
    }
  }
}

warn("keys are ", join(",", keys %topu1), "\n");

#find matches between top and bottom strand
for my $ctg ( sort keys %topu1 ) {
  warn("there are ", (scalar keys %{$topu1{$ctg}}), " locations in $ctg\n");
  for my $location (sort { $a <=> $b } keys %{$topu1{$ctg}} ) {
    my $location2 = $location + 9;
    if (exists $botu1{$ctg}->{$location2}) {
      chop $topu1{$location};
      chop $botu1{$location2};
      print join("\t", ($location, $location2,$topu1{$location},$botu1{$location2})),"\n";
    }
  }
}
