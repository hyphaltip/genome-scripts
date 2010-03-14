#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::Seq::Quality;

my $input_base = shift;
my $base = $input_base;
$base =~ s/\.(fasta|fa)$//;

my $fain = Bio::SeqIO->new(-format => 'fasta',
			   -file   => $input_base);
my $qualin= Bio::SeqIO->new(-format => 'qual',
			    -file   => $input_base .".qual");

my $fastqout = Bio::SeqIO->new(-format => 'fastq-solexa',
			       -file   => ">$base.fq");
			    
# $data is a hash reference containing all arguments to be passed to
# the Bio::Seq::Quality constructor
while (my $seq = $fain->next_seq) {
    my $qual = $qualin->next_seq;
    if( $qual->id ne $seq->id ) {
	die("IDs do not match, sequences are not in same order in $input_base and $input_base.qual\n");
    }
    # process $data, such as trim, etc
    my $data = Bio::Seq::Quality->new(-qual => $qual,
				     -id   => $qual->id,
				     -seq  => $seq->seq,
				     );
    $fastqout->write_seq($data);
}

