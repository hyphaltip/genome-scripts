#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $head = <>;
my ($chrom,$pos,$ref,@strains) = split(/\s+/,$head);
my $total = scalar @strains;
my %seqs;
while(<>) {
    my ($chr,$p,$ref_allele, @alleles) = split;
    my %count = ( $ref_allele => 1 );
    my @sanitized;
    for my $allele ( @alleles ) {
	my ($al1,$al2) = split(/\//,$allele);
	if( $al2 ) {
	    # this is diploid
	    # do something different?
	    die " cannot process diploid tab files for now";
	} 
	$al1 =~ s/\./-/;
	$count{$al1}++;
	push @sanitized, $al1;
    }

    if( exists $count{'-'} && 
	(my $per = $count{'-'} / $total) > 0.8 ) {
	warn(sprintf("skipping $chr,$p as %.2f%% are uncalled\n", 
		     $per));
	next;
    }
    $seqs{'REF'} .= $ref_allele;
    my $n = 0;
    for my $al1 ( @sanitized ) {
	$seqs{$strains[$n++]} .= $al1;	
    }
}
my $out = Bio::SeqIO->new(-format => 'fasta',-file => ">vcfseqs.fas");
for my $strain ( keys %seqs ) {
    $out->write_seq(Bio::Seq->new(-id => $strain,
				  -seq => $seqs{$strain}));
}
