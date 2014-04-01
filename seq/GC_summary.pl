#!/usr/bin/perl -w
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(sum);

my $format = 'fasta';

GetOptions('f|format:s' => \$format);

my $infile = shift;
my $in = Bio::SeqIO->new(-format => $format, 
			 -file => $infile);

my $cumul = 0;
my $len;
while(my $s = $in->next_seq ) {
    my $gc = ($s->seq =~ tr/GCgc/GCgc/);
    $cumul += $gc;
    $len += $s->length;
}

printf "%.2f%% GC \n",100 * $cumul/$len;
