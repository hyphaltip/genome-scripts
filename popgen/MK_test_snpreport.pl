#!/usr/bin/perl -w
# $Id: MK_test_snpreport.pl 249 2008-12-12 22:20:01Z stajich $

# this is intended specifically for Coccioidies project
# processing the output from Broad SNP calling pipeline

use strict;
use Getopt::Long;
use List::Util qw(sum);
use Text::NSP::Measures::2D::Fisher2::right;
use Text::NSP::Measures::2D::Fisher2::twotailed;
use Text::NSP::Measures::2D::Fisher2::left;

my $min_AF = 0.25;
my $reportfile = 'SnpReport_Cspp-3.txt';
my $refallele_header = 'a_CI2';
my ($outerr,$outfile);
GetOptions(
	   'o|out:s' => \$outfile,
	   'e|err:s' => \$outerr,
	   'i|in:s'  => \$reportfile,
	   'maf=f'   => \$min_AF,  
	   );

# open input filehandle
open(my $fh => $reportfile) || die $!;
$outerr ||= "$reportfile.MK_ambiguous.dat";
$outfile ||= "$reportfile.MK.dat";
open(my $errfh => ">$outerr") || die $!;
open(my $outfh => ">$outfile") || die $!;

my $header = <$fh>; # first line is header, read it in
 
chomp($header); # drop last character
my @header = split(/\t/,$header); # turn it into a list
# turn list into hash where the keys are the col names
# and values are the col number
my $i =0;
my (%headerhash,@cicols,@cpcols);
for my $col ( @header ) {
    # if col starts with a_CI it is Cimmitis, a_CP it is cposadasii
    if( $col =~ /^a_CI/ ) {
	push @cicols, $i;
    } elsif( $col =~ /^a_CP/ ) {
	push @cpcols, $i;
    }
    $headerhash{$col} = $i++;
}
my %genes;
while(<$fh>) {
    chomp; # drop newline
    my @row = split(/\t/,$_); # row out of a line
    my ($snp) = $row[0];
    # want only the rows which aren't intergenic so skip these
    next unless ( $row[ $headerhash{'snptype'} ] eq 'Synonymous' ||
		  $row[ $headerhash{'snptype'} ] eq 'Non-Synonymous' );
#     want only those would good quality
#    next if ( $row[ $headerhash{'translation_quality'} ] ne 'good' );
    my $genename = $row[ $headerhash{'genes'} ];
    if( $genename =~ /\,\s/) {
#	warn("gene has comma or space ($genename)!");
	next;
    }
    my @bits; # bit one is within (0) or between (1)
              # bit zero  is syn    (0) or nonsyn  (1)
    # 2x2 table will be 
    # 0 00 -> Syn    within  (polymorphism)
    # 1 01 -> NonSyn within  (polymorphism)
    # 2 10 -> Syn    between (fixed)
    # 3 11 -> NonSyn between (fixed)
    
    $bits[0] = $row[ $headerhash{'snptype'} ] eq 'Non-Synonymous' ? 1 : 0;

    my $refallele = $row[ $headerhash{$refallele_header} ];
    my (%ci_alleles,%all_alleles);
    for my $ci ( @cicols ) {
	next if $row[$ci] eq '' || length($row[$ci]) == 0;
	# this is going to count the number of each type of allele
	$ci_alleles{$row[$ci]}++;
	# keep track of total alleles for all individuals
	$all_alleles{$row[$ci]}++;
    }
    my %cp_alleles;
    for my $cp ( @cpcols ) {
	next if $row[$cp] eq '' || length($row[$cp]) == 0;
	# this is going to count the number of each type of allele
	$cp_alleles{$row[$cp]}++;
	# keep track of total alleles for all individuals
	$all_alleles{$row[$cp]}++;
    }
    my $total_alleles = scalar keys %all_alleles;
    if( $total_alleles > 2 ) {
	print $errfh join("\t", $snp, $genename,"TOO MANY ALLELES", join(",", sort keys %all_alleles)),"\n";
	next;
    }
    my @cialleles = sort keys %ci_alleles;
    my $ci_allele_count = scalar @cialleles;
    my @cpalleles = sort keys %cp_alleles;
    my $cp_allele_count = scalar @cpalleles;
    
    # this would be a good point to skip SNPs where the freq is < 0.25 
    # or something, use the mAF column
    # next if $row[ $headerhash{'mAF'} ] < $min_AF;
    
    if( $cp_allele_count == 1 && $ci_allele_count == 1 ) {	
	# monomorphic for both so this is a fixed different
	$bits[1] = 1;
    } else {
	$bits[1] = 0;	
    }    
    my $index = $bits[0] + 2*$bits[1];
    $genes{ $genename }->[$index]++;
}
print $outfh join("\t", qw(GENE SYN_POLYMORPH NONSYN_POLYMORPH 
			   SYN_FIXED NONSYN_FIXED 
			    PVALUE_LEFT PVALUE_RIGHT PVALUE_TWOTAILED)),"\n";

for my $gene ( sort keys %genes ) {
    for my $i ( 0..3 ) {
	$genes{$gene}->[$i] ||=0;
    }
    
# according to Text::NSP::Measures::2D::Fisher::twotailed
# this is how we format the data for the 2x2 table for input
# into calculateStatistic
#                  fixed     poly
#        N         n11      n12 ? n1p
#        S         n21      n22 ? n2p
#                  --------------
#                  np1      np2   npp

    my $pvalue_right = &Text::NSP::Measures::2D::Fisher2::right::calculateStatistic
	(n11=> $genes{$gene}->[3], # N_fixed
	 n1p=> $genes{$gene}->[3] + $genes{$gene}->[1], # N_fixed + N_poly
	 np1=> $genes{$gene}->[3] + $genes{$gene}->[2], # N_fixed + S_fixed
	 npp=> sum(@{$genes{$gene}})  # N_poly + N_fixed + S_poly + S_fix
	 );	 

    my $pvalue_left = &Text::NSP::Measures::2D::Fisher2::left::calculateStatistic
	(n11=> $genes{$gene}->[3], # N_fixed
	 n1p=> $genes{$gene}->[3] + $genes{$gene}->[1], # N_fixed + N_poly
	 np1=> $genes{$gene}->[3] + $genes{$gene}->[2], # N_fixed + S_fixed
	 npp=> sum(@{$genes{$gene}})  # N_poly + N_fixed + S_poly + S_fix
	 );	 

    my $pvalue_twotail = &Text::NSP::Measures::2D::Fisher2::twotailed::calculateStatistic
	(n11=> $genes{$gene}->[3], # N_fixed
	 n1p=> $genes{$gene}->[3] + $genes{$gene}->[1], # N_fixed + N_poly
	 np1=> $genes{$gene}->[3] + $genes{$gene}->[2], # N_fixed + S_fixed
	 npp=> sum(@{$genes{$gene}})  # N_poly + N_fixed + S_poly + S_fix
	 );	 

    print $outfh join("\t",$gene,@{$genes{$gene}}, 
		      sprintf("%.3g",$pvalue_left),
		      sprintf("%.3g",$pvalue_right),
		      sprintf("%.3g",$pvalue_twotail),
		      ),"\n";

}
