#!/usr/bin/perl -w
use strict;

use Bio::DB::Sam;
use Getopt::Long;

my $debug = 0;
my $features;
GetOptions(
	   'v|verbose!' => \$debug,
	   'f|gff|feature:s' => \$features,
	   );
my ($bam,$fa) = @ARGV;
# high level API
my $sam = Bio::DB::Sam->new(-bam  =>$bam,
			    -fasta=>$fa,
			    );
my @bases = qw(C A G T);

my %skip;
if( $features ) {
    open(my $fh => $features ) || die $!;
    while(<$fh>) {
	next if /^\#/;
	my ($seqid,$src,$type,$start,$end) = split;
	push @{$skip{$seqid}}, [ $start, $end];
    }
}
my %counts;
my @targets    = $sam->seq_ids;
for my $target ( @targets ) {
    my $iterator   = $sam->features(-iterator=>1, -seq_id => $target);

    while( my $aln = $iterator->next_seq ) {
	my $seqid  = $aln->seq_id;
	my $start  = $aln->start;
	my $end    = $aln->end;
	#my $strand = $aln->strand;
	#my $cigar  = $aln->cigar_str;
	my $len = $aln->length;
	warn join("\t", $seqid, $aln->name, $len),"\n" if $debug;
	my $dna = $aln->query->dna ;
	my %f;
	for ( split('',$dna) ) {
	    $f{$_}++;
	}
	next if keys %f <= 2; # drop those AAA or TTT runs
	my $skip_this = 0;
	if( exists $skip{$seqid} ) {
	    for my $l ( @{$skip{$seqid}} ) {
		if( &overlaps($start,$end,@$l) ) {
		    $skip_this = 1;
		    last;
		}
	    }
	}
	next if $skip_this;
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

sub overlaps {
    my ($f1_s,$f1_e,
	$f2_s,$f2_e) = @_;
    return not ( $f1_s > $f2_e ||
		 $f1_e < $f2_s);
}
