#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
my $merged_indels = shift;

my $sep = ';';
open(my $fh => $merged_indels) || die $!;
my %strains;
my $odir = 'perchrom';
mkdir($odir) unless -d $odir;
my %data;
while(<$fh>) {
    chomp;
    
    my ($orf,$strain,$contig,$pos,$size_seq,$expr) = split($sep,$_);
    $data{$contig}->{$pos}->{$strain} = $size_seq;
    $strains{$strain}++;
}
my @strains = sort { $a <=> $b} keys %strains;
print join(",", @strains),"\n";
my %all;
for my $chrom ( keys %data ) {

#    my $alnio = Bio::AlignIO->new(-format => 'nexus',
#				  -file   => ">$odir/$chrom.indels.nex");
    open(my $ofh => ">$odir/$chrom.indels.phy") || die $!;
    my %d;
    my $count = 0;
    for my $pos ( sort { $a <=> $b } keys %{$data{$chrom}} ) {
	my %seen;
	for my $s ( @strains ) {
	    my $allele = exists $data{$chrom}->{$pos}->{$s} ? $data{$chrom}->{$pos}->{$s} : '0';
	    push @{$seen{$allele}}, $s;
	}
	# check the data before adding it
	if( keys %seen > 2 ) { # more than 1 insert size seen
	    warn("multiple insertion size alleles at $chrom $pos\n");
	    next;
	}
	while( my ($allele,$strains) = each %seen ) {
	    for my $s ( @$strains ) {
		$d{$s} .= $allele ? 1 : 0;
	    }
	}
	$count++;
    }
    
 #   my $aln = Bio::SimpleAlign->new;
    printf $ofh "%7d %d\n",scalar @strains, $count; 
    for my $strain ( @strains ) {
	$all{$strain} .= $d{$strain};
	printf $ofh "str%-10s %s\n",$strain,$d{$strain};
    }

#    $alnio->write_aln($aln);
}

open(my $allfh => ">all.indels.phy") || die $!;
printf $allfh "%7d %d\n",scalar @strains, length $all{$strains[0]};
for my $strain ( @strains ) {    
    printf $allfh "str%-10s %s\n",$strain,$all{$strain};
}

