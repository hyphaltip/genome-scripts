#!/usr/bin/perl -w
use strict;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;

my $seqfile;
GetOptions('s|seqfile:s' => \$seqfile);

$seqfile ||= shift @ARGV;
my $db = Bio::DB::Fasta->new($seqfile);

my $out = Bio::SeqIO->new(-format => 'fasta');
while(<>) {
    next if /^\#/;
    my ($query,$hit,$pid,$aln_len, $mm, $gap_open,
	$qstart, $qend, $sstart, $send, 
	$evalue, $bitscore) = split;
    if( my $seq = $db->seq($hit,$sstart,$send) ) {
	$out->write_seq(Bio::Seq->new(-id => sprintf("%s_%s-%d",
						     $hit,$sstart,$send),
				      -seq => $seq,
				      -description => 
				      sprintf('PID=%s AlnLen=%d Evalue=%s Bits=%s QueryAln=%d..%d',
					      $pid,$aln_len,$evalue,$bitscore,
					      $qstart,$qend)));
    }
}
