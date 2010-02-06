#!/usr/bin/perl -w
use strict;

use Bio::DB::Sam;
use Getopt::Long;
my $debug = 0;
GetOptions(
	   'v|verbose!' => \$debug,
	   );
my ($bam,$fa) = @ARGV;
# high level API
my $sam = Bio::DB::Sam->new(-bam  =>$bam,
			    -fasta=>$fa,
			    );
my @bases = qw(C A G T);
my %counts;
my @targets    = $sam->seq_ids;
for my $target ( @targets ) {
    my $iterator   = $sam->features(-iterator=>1, -seq_id => $target);

    while( my $aln = $iterator->next_seq ) {
	my $seqid  = $aln->seq_id;
	#my $start  = $aln->start;
	#my $end    = $aln->end;
	#my $strand = $aln->strand;
	#my $cigar  = $aln->cigar_str;
	my $len = $aln->length;
	warn join("\t", $seqid, $aln->name, $len),"\n" if $debug;
	my $dna = $aln->query->dna ;
	my %f ;
	for ( split('',$dna) ) {
	    $f{$_}++;
	}
	next if keys %f <= 2;
	my $base = substr($dna,0,1);
	if( $aln->name =~ /^(\d+)\-(\d+)$/) {
	    $counts{$len}->{$base} += $2;
	} else {
	    $counts{$len}->{$base}++;
	}
	last if $debug;
    }
    last if $debug;
}
print join("\t",qw(LENGTH), @bases), "\n"; 
for my $size ( sort { $a <=> $b } keys %counts ) {    
    print join("\t", $size, map { $counts{$size}->{$_} || 0 } @bases),"\n";
}
