#!/usr/bin/perl -w
use strict;

my $file = shift || die "cannot run without input file";

open(my $fh => $file) || die "cannot open $file\n";
my %types;
while (<$fh>) {
  next if /^\#/ || /^\s+$/;;
  my ($chrom,$src,$type,$start,$end,$score,$strand) = split;
  for my $pos ( $start..$end ) {
    $types{$chrom}->{$pos}->{$type}++;
  }
}

my %totals;
for my $chrom ( keys %types ) {
  for my $t ( values %{$types{$chrom}} ) {
    for my $kind ( keys %$t ) {
# kind = gene, CDS, tRNA, etc
      $totals{$kind}++;
    }
  }
}
for my $type ( keys %totals ) {
  print "$type $totals{$type} nt\n";
}
