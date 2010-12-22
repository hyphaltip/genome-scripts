#!/usr/bin/perl -w
use strict;

use Bio::DB::Sam;

my $debug = 0;

my ($base,$genome,$regions,$odir) = @ARGV;
$odir ||= "reads_per_region";
mkdir($odir);
my $db1 = Bio::DB::Sam->new(-bam => "$base.1.bam",
			   -fasta => $genome);

my $db2 = Bio::DB::Sam->new(-bam => "$base.2.bam",
			   -fasta => $genome);

open(my $fh => $regions ) || die $!;
unlink("$odir/$base.1.hits");
unlink("$odir/$base.2.hits");
while(<$fh>) {
    my ($chr,$source,$type,$start,$end,undef,$strand,undef,$name) = split;
    if( $name =~ /(ID|Name)=(\S+)/i ) {
	$name = $2;
    }

    my $segment = $db1->segment($chr,$start => $end);
    my $iterator = $segment->features(-iterator => 1);
    open(my $ofh => ">$odir/$base.$name.reads1") || die $!;
    while( my $aln = $iterator->next_seq ) {
	my $rname = $aln->query->seq_id;
	my $hit = `grep $rname $base.2.sam | tee -a $odir/$base.2.hits`;
	my @other_nfo;
	if( $hit ) {
	    chomp($hit);
	    warn("hit is $hit\n") if $debug;
	    my ($id,$flag,$chrom,$start,$length,$cigar) = split(/\s+/,$hit);
	    if( $chrom ne '*' ) {
		my $end = $start + $length;
		if( $flag & 0x10) {
		    $end = $start - $length;
		} 
		@other_nfo = ($chrom,$start,$end,$cigar);
	    }
	}
	print $ofh join("\t", 
			$aln->seq_id,
			$aln->start,
			$aln->end,			
#			$aln->dna,
			$aln->query->seq_id,
			#$aln->query->start,
			#$aln->query->end,
			#$aln->query->strand,			
			$aln->query->dna,
			$aln->cigar_str,
			@other_nfo
			),"\n";
    }

    $segment = $db2->segment($chr,$start => $end);
    $iterator = $segment->features(-iterator => 1);
    open($ofh => ">$odir/$base.$name.reads2") || die $!;
    while( my $aln = $iterator->next_seq ) {
	my $rname = $aln->query->seq_id;
	my $hit = `grep $rname $base.1.sam | tee -a $odir/$base.1.hits`;
	my @other_nfo;
	if( $hit ) {
	    chomp($hit);
	    my ($id,$flag,$chrom,$start,$length,$cigar) = split(/\s+/,$hit);
	    if( $chrom ne '*' ) {
		my $end = $start + $length;
		if( $flag & 0x10) {
		    $end = $start - $length;
		} 
		@other_nfo = ($chrom,$start,$start+$length,$cigar);
	    }
	    warn("hit is $hit\n") if $debug;
	}
	print $ofh join("\t", 
			$aln->seq_id,
			$aln->start,
			$aln->end,			
#			$aln->dna,
			$aln->query->seq_id,
#			$aln->query->start,
#			$aln->query->end,
#			$aln->query->strand,
			
			$aln->query->dna,
			$aln->cigar_str,
			@other_nfo,
			),"\n";
    }
}
