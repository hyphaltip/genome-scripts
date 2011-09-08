#!/usr/bin/perl -w
use strict;
my $file = shift;
my $libsize = shift;
open(my $fh => "grep '^>' $file |") || die $!;

while(<$fh>) {
 if(/^>(\S+)/) {
  my $id = $1;
  if( $id =~ /(\S+)\/[1]/ ) { 
   my $nm = $1;
   print "$nm/1 $nm/2 $libsize\n"; 
 } else {
  warn("unable to process, expect only 1 kind of read naming now ($id)\n");
 }
}
}
