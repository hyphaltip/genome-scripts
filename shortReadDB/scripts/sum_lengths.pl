#!/usr/bin/perl -w
use strict;
my %types;
while (<>) {
  next if /^\#/ || /^\s+$/;
  my ($seqid,$src,$type,$start,$end) = split;
  $types{$type} += abs($end - $start) + 1;
}

print join("\t",qw(TYPE COUNT)),"\n";
for my $type ( sort keys %types) {
  print join("\t", $type,$types{$type}),"\n";
}
