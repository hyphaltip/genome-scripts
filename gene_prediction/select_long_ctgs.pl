#!env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $min_len = 100_000;

GetOptions( 'm|min:s' => \$min_len );
my $in = Bio::SeqIO->new(-format => 'fasta',
			 -file   => shift);

my $out = Bio::SeqIO->new(-format => 'fasta');
while ( my $s = $in->next_seq ) {
    next unless $s->length > $min_len;
    $out->write_seq($s);
}
