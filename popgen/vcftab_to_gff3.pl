#!env perl
use strict;
use warnings;

my $in = shift || die $!;

open(my $fh => $in) || die $!;
my @header;
while(<$fh>) {
    if(s/^#CHROM/CHROM/) {
	@header = split;	
    } else {
	my ($chrom,$pos,$ref,@strains) = split;
	my @grp;
	print join("\t", $chrom,qw(GATK SNP),$pos, $pos, '.','+','.',sprintf("ID=%s;alleles=%s;ref_allele=%s;");
    }
}
