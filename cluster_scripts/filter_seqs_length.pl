#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => shift);
my $out = Bio::SeqIO->new(-format => 'fasta');
my $min_size = shift || 10_000;
while( my $seq = $in->next_seq ) {
  $out->write_seq($seq) if $seq->length >= $min_size; 
}
