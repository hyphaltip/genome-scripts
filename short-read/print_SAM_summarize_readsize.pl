#!/usr/bin/perl -w
use strict;

use Bio::DB::Sam;

my $bamfile = shift;
my $name = $bamfile;
$name =~ s/\.bam$//;

my $dbh = Bio::DB::Sam->new(-bam => $bamfile);

# get all the reads, grab 5' base, make frequency plot

my @bases = qw(C A G T);

for my $seq ( map { $db->segment(-seq_id => $_) } @targets ) {
	print $seq->seq_id, " ", $seq->length, " ", $seq->start, "..",$seq->end,"\n";
        for my $read ( $seq->features ) {
	    my $length = $read->query->length;
	    my $ofh = $ofh{$length};
	    unless( defined $ofh ) {
		open($ofh => ">$dir/$base.$length.dat") || die $!;
		$ofh{$length} = $ofh;
	    }
for my $target ( $sam->seq_ids ) {
    for my $aln ( $dbh->get_features_by_location(-seq_id => $target) ) {
	my $qstart  = $aln->start;
	my $qend    = $aln->end;
	my $qstrand = $aln->strand;
	my $len = $aln->length;
	my $dna = $aln->query->dna;
	my %f;
	for ( split('',$dna) ) {
	    $f{$_}++;
	}
	next if keys %f <= 2; # drop those AAA or TTT runs
	my $skip_this = 0;
	
    }
}
