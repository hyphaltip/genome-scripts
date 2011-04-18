#!/usr/bin/perl -w
use strict;

# input file is .singletons file from vcftools --singletons option
my %strains;
while(<>) {
    next if /INDV/;
    my ($chrom,$pos,$type,$allele,$strain) = split;
    $strains{$strain}->{$type}++;
}
print join("\t", qw(STRAIN DOUBLES SINGLES)),"\n";
for my $strain ( keys %strains ) {
    print join("\t", $strain, $strains{$strain}->{D} || 0, $strains{$strain}->{S} || 0),
    "\n"; 
}
