#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw(shuffle);
use Getopt::Long;

my @header;
my %inds;
my @loci;

my $missing_data = -9;
my $rand_sample = 5000;
my $sample2popfile;

warn("USE plink instead and faststrutucture. This should be deprecated");

GetOptions(
    'p|pop:s'   => \$sample2popfile,
    'rs|rand:i' => \$rand_sample,
    );

while(<>){
    if(/^#?CHROM/ ){
	s/^#//;
	@header = split;
    } elsif( /^#/ ) { next } 
    else {
	my ($chrom,$pos,@gt) = split;
	push @loci, sprintf("%s.%d",$chrom,$pos);
	my $i =2 ; # skip first 2 cols of 
	for my $gt ( @gt ) {
	    $gt = $missing_data if $gt eq 'NA';
	    my ($gtA,$gtB) = ($gt,$gt);
	    if( $gt eq '1' ) {
		$gtA = 0;
		$gtB = 2;
	    }
	    $inds{$header[$i++]}->{$loci[-1]} = [$gtA,$gtB];
	}
    }
}
my %sample2pop;
if( $sample2popfile ) {
    open(my $fh => $sample2popfile) || die "$sample2popfile: $!";
    while(<$fh>) {
	my ($sample,$pop) = split;
	$sample2pop{$sample} = $pop;
    }
}

my ($chrom,$pos,@names) = @header;
my @keep_loci = &randsample($rand_sample,\@loci);
print join(" ", @keep_loci), "\n";

for my $ind ( @names ) {    
    my $pop = $sample2pop{$ind} || 0;    
# FIX ME HERE
    print join(" ", $ind, $pop, map { @{$inds{$ind}->{$_}} } @keep_loci),"\n";    
    #print join(" ", $ind, $pop, map { $inds{$ind}->{$_} } @keep_loci),"\n";    
}

sub randsample {
    my ($count,$loci) = @_;
    my @v = shuffle @$loci;
    splice(@v,0,$count);
}
