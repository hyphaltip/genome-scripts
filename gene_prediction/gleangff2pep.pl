#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $out = Bio::SeqIO->new(-format => 'fasta');
while(<>) {
    chomp;
    next if /^\#/ || /^\s+$/);
    my @line = split(/\t/,$_);
    if( $line[2] eq 'mRNA' && 
	$line[-1] =~ /ID=Gene:([^;]+);translationSeq=([^;]+);/ ) {
	$out->write_seq(Bio::Seq->new(-id => $1,
				      -seq => $2));
    }
}
