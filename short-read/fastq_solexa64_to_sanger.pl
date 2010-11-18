#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Seq::Quality;

my $in = Bio::SeqIO->new(-format => 'fastq-solexa',
			 -fh     => \*ARGV);
my $out = Bio::SeqIO->new(-format => 'fastq');


# $data is a hash reference containing all arguments to be passed to
# the Bio::Seq::Quality constructor
while (my $data = $in->next_dataset) {
    $out->write_seq(Bio::Seq::Quality->new(%$data));
}

