#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
use Statistics::Descriptive;

# this script expects a GFF file as input
# and that it is GFF3 with Parent and ID tags

my $file = shift || die "cannot run without input file";

open(my $fh => $file) || die "cannot open $file\n";
my %types;
while (<$fh>) {
  next if /^\#/ || /^\s+$/;;
  my ($chrom,$src,$type,$start,$end,$score,$strand,$frame,$info) = split;
  my %info = map { split(/=/,$_) } split(/;/,$info);
  my ($parent) = ($info{'Parent'} || $info{'ID'} );
  $types{$type}->{$parent} += abs($end - $start);
}

for my $type ( keys %types ) {
  my @lens;
  for my $gene ( keys %{$types{$type}} ) {
   push @lens, $types{$type}->{$gene};
  }
 my $stats = Statistics::Descriptive::Full->new();
 $stats->add_data(@lens);
 printf "%s mean=%.2f median=%d total_count=%d total_length=%d\n",$type, $stats->mean, $stats->median, $stats->count,$stats->sum;
}
