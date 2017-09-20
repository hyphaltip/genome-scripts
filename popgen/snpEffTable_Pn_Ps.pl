#!/usr/bin/perl
use strict;
use warnings;

my $INFILE = 'snpEff_genes.txt';
if( @ARGV ) {
    $INFILE = shift @ARGV;
}

print join("\t",qw(Gene Pn Ps Pn_Ps Total_Variants)),"\n";
open(my $IN => $INFILE) || die "Cannot open $INFILE";

my $header;
while( <$IN> ) {
    if( /^\#GeneName/ ) {
	$header = $_;
	last;
    }
}
if ( ! $header ) {
    die("Header is not in correct format, expected a starting '#'\n");
}
# PARSE HEADER
my $i =0;
my %cols;

%cols = map { $_ => $i++ } split(/\s+/,$header);
#warn("col names are:\n", join("\n",sort { $cols{$a} <=> $cols{$b} } keys %cols), "\n");

# PARSE FILE
while (<$IN>) { 
    my @row = split(/\s+/,$_);
    # the columns names are keys in the %cols hash
    # this returns the index in the array @row
    # which is just split on whitespace
    # this way we don't have to remember the order of the columns
    # and this is readable
    my $Pn = $row[$cols{variants_effect_missense_variant}] + 
	     $row[$cols{variants_effect_start_lost}] + 
	     $row[$cols{variants_effect_stop_gained}] + 
	     $row[$cols{variants_effect_stop_lost}];

    my $Ps = $row[$cols{variants_effect_synonymous_variant}] + 
 	     $row[$cols{variants_effect_stop_retained_variant}];
    my $Total=$Pn + $Ps;
     my $PnPsRatio = ($Ps == 0 ) ? 'Undef' : sprintf("%.4f",$Pn/$Ps); # avoid divide by 0 issues
    print join("\t",$row[$cols{GeneId}],
	       $Pn,$Ps,$PnPsRatio,$Total), "\n";
}
