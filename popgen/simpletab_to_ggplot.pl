#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $binsize = 10_000;
my $odir = 'freq-per-chrom';
GetOptions(
	   'b|bin:i' => \$binsize,
	   'o|odir:s'=> \$odir,
	   );
mkdir($odir) unless -d $odir;
chomp(my $hdr = <>);
my ($CHR,$POS,$REF,@strains) = split(/\t/,$hdr);

my %bins;
my $bincount = 0;

while(<>) {
    my ($chrom,$pos,$ref,@row) = split;
    my $bin = int($pos / $binsize);
    for my $strain ( @strains ) {
	my $allele = shift @row;
	my @alleles = split '',$allele;
	my %count;
	for my $al ( @alleles) { $count{$al}++ }
	my $type = keys %count == 1 ? 'hom' : 'het';
	next if( $type eq 'hom' && $count{$ref} );
	$bins{$chrom}->{$strain}->{$bin}->{$type}++;
    }
}

for my $chrom ( keys %bins ) {
    open(my $fh => ">$odir/SC\_$chrom.dat") || die $!;
    print $fh join("\t", qw(BIN SNPCOUNT TYPE STRAIN)), "\n";
    for my $strain ( sort  keys %{$bins{$chrom}} ) {
	for my $bin ( sort { $a <=> $b } keys %{$bins{$chrom}->{$strain}} ) {
	    for my $type ( qw(hom het) ) {
		print $fh join("\t", $bin * $binsize, 
			       $bins{$chrom}->{$strain}->{$bin}->{$type} || 0,
			       uc($type), $strain), "\n";
	    }
	}
    }
}
			     
