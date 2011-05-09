#!/usr/bin/perl -w
use strict;

my %types;
while (<>) {
  next if /^\#/;
  my ($chrom,$src,$type,$start,$end,$score,$strand) = split;
  for my $pos ( $start..$end ) {
    $types{$type}->{$chrom}->{$pos}++;
  }
}

for my $type ( keys %types ) {
  my $total = 0;
  for my $chr ( keys %{$types{$type}} ) {
    # count the number nt that are covered by this feature type;
    $total += scalar keys %{$types{$type}->{$chr}};
  }
  print "$type $total nt\n";
}
