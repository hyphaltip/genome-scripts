#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $species = "Coprinus_cinereus";
GetOptions(
	's|species:s'=> \$species);
mkdir($species);
mkdir("$species/data");
my @lines; 
while(<>){ 
next if /^\#/;
 my @line = split; 
 if ( $line[2] eq "mRNA" ) { 
 $line[-1] =~ s/ID=([^;]+);\S+/$1/; 
 push @lines, [$line[0], $line[3], $line[-1]];
 } 
} 

my $last;
my $fh;
open(my $name => ">$species/data_list") || die $!;
for my $l ( sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @lines ) { 
 if( ! defined $last || $l->[0] ne $last ) {
   close($fh) if defined $fh;
   print $name "$species/data/".$l->[0].".data", "\n";
   open($fh, ">$species/data/".$l->[0].".data") || die $!;
 }   
 print $fh join("\t", $l->[2], $l->[1], $l->[0]), "\n"; 
 $last = $l->[0];
}
close($name);
close($fh);
