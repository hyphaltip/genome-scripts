#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;
use Getopt::Long;
my $format = 'fasta';
my $length = 100;
GetOptions("f|fasta:s" => \$format,
	   'l|len|length:i' => \$length);

my $in = Bio::SeqIO->new(-fh => \*ARGV,
			 -format=>$format);
my $out = Bio::SeqIO->new(-format => $format);
while(my $s = $in->next_seq ){
    next unless $s->length >= $length;
    $out->write_seq($s);
}
