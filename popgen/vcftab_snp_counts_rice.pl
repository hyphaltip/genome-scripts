#!/usr/bin/perl -w
use strict;
use warnings;
my $header = <>;
my (undef,undef,undef,@strains) = split(/\s+/,$header);
my %counts;
while(<>) {
    my ($chr,$pos,$ref,@alleles) = split;
    my $i = 0;
    my (%strain,%all_alleles);
    for my $genotype ( @alleles ) {
	my @al = split(/\//,$genotype);
	my $strain = $strains[$i++];
	for ( @al ) {
	    $strain{$strain}->{$_}++;
	    $all_alleles{$_}++;
	}
    }
    next if( exists $all_alleles{'.'} );# should we skip the whole position?
    for my $st ( @strains ) {
	{}
    if( ! exists $all_alleles{$ref} ) {
	$counts{'shared-derived'}++;
    }
}
