#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $min_diff = 4;
my $orig = Bio::SeqIO->new(-format => 'fasta',
			   -file   => shift @ARGV);


my $trim = Bio::SeqIO->new(-format => 'fasta',
			   -file   => shift @ARGV);

my %len;
while(my $seq = $orig->next_seq ) {
    $len{$seq->id} = $seq->length;
}

while(my $seq = $trim->next_seq ) {
    if( abs($seq->length - $len{$seq->id}) > $min_diff ) {
	print join("\t", $seq->id, $len{$seq->id}, $seq->length), "\n";
 }   
}

