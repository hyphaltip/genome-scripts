#!/usr/bin/perl -w
use strict;
use warnings;
my $header = <>;
my (undef,undef,undef,@strains) = split(/\s+/,$header);
my %counts;
my %hets;
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
	if( keys %{$strain{$strain}} == 2 ) {
	    $hets{$strain}++;
	}
    }

    next if( exists $all_alleles{'.'} ); # should we skip the whole position?
    
    for my $st ( @strains ) {
	if( ! defined $strain{$st}->{$ref} ) {
	    # doesn't have the NB ref allele	    
	    $counts{'diff_ref-'.$st}++;# += 2; # since neither allele matches add 2
	} else {
	    $counts{'shared_ref-'.$st}++;# += $strain{$st}->{$ref};
	}	
    }
    for( my $j =0; $j < @strains; $j++ ) {
	my $straina = $strains[$j];
	for( my $k = $j+1; $k < @strains; $k++ ) {
	    my $strainb = $strains[$k];

	    for my $allele_A (  keys %{$strain{$straina}} ) {
		if( ! exists $strain{$strainb}->{$allele_A} ) {
		    $counts{"diff_$straina-$strainb"}++;# += $strain{$straina}->{$allele_A};
		} else {
		    $counts{"shared_$straina-$strainb"}++;# += $strain{$straina}->{$allele_A};
		}
	    }
	    for my $allele_B ( keys %{$strain{$strainb}} ) {
		if( ! exists $strain{$straina}->{$allele_B} ) {
		    # keep strains in this order so it will be simpler
		    $counts{"diff_$straina-$strainb"}++;# += $strain{$strainb}->{$allele_B};
		} 
		# already counted similarities, would have been in the previous count
	    }
	}
    }
}

for my $strain (sort keys %hets ){
    print join("\t", 'HET', $strain,$hets{$strain}),"\n";
}

for my $kind ( grep { /diff/ } sort keys %counts ) {
    print join("\t", 'DIFF', $kind, '',$counts{$kind} ),"\n";
}

for my $kind ( grep { /shared/ } sort keys %counts ) {
    print join("\t", 'SHARED', $kind,'', $counts{$kind} ),"\n";
}

# what about arabidposis het positions?

