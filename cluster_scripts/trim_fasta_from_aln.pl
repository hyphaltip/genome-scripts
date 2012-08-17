#!/usr/bin/perl
use warnings;
use strict;

use Bio::SeqIO;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => shift);
my ($qid,$qlen,$seq,$ct);
my %match;
while(<>) {
    if(/^\#\s+Query:\s+(\S+)\s+-\s+(\d+)\s+nt/ ) {
	($qid,$qlen) = ($1,$2);	
	$seq = $in->next_seq;
	if( $seq->id ne $qid ) {
	    warn("db and report are not in sync: $qid != ", $seq->id,"\n");
	    last;
	}
	$ct = 0;
    } elsif(/^\#/) { next }
    else {
	my ($query,$hit,$pid,
	    $identbases,$mm,$gap,
	    $qstart,$qend,$hstart,$hend,$hevalue,$hbit) = split;
	my $strand = ( $qstart < $qend ) ? 1 : -1;
	if( $query ne $qid ) {
	    warn("query id ($query) is not same as $qid\n");
	    last;
	}
	if( $query ne $seq->id ) {
	    warn("query id ($query) is not same as seq id ", $seq->id,"\n");
	}
	($qstart,$qend) = sort { $a <=> $b } ($qstart,$qend);
	
	my ($clip_left, $clip_right) = ([1,1],[$seq->length,$seq->length]);
	if( $qstart != 1 ) {
	    $clip_left = [1,$qstart];
	}
	if( $qend < $seq->length ) {
	    $clip_right = [$qend,$seq->length];
	}
	my $keep = [$clip_left->[1],
		     $clip_right->[0]];

	if( $ct++ == 0 ) {
	    if( $strand > 0 ) {
		print join("\t", $query, $strand, join("-", @$clip_left),
			   join("-", @$clip_right),join("-",@$keep),
			   $clip_left->[1] > 1 ? $seq->subseq(@$clip_left) : '.',		       
			   ($clip_right->[0] < $seq->length) ? $seq->subseq(@$clip_right) : '.',
			   $seq->subseq(@$keep),
		    ),"\n";
	    } else {
		print join("\t", $query, $strand, 
			   join("-",reverse @$clip_left),
			   join("-",reverse @$clip_right),
			   join("-",reverse @$keep),
			   $clip_left->[1] > 1 ? $seq->trunc(@$clip_left)->revcom->seq : '.',
			   ($clip_right->[0] < $seq->length) ? $seq->trunc(@$clip_right)->revcom->seq : '.',
			   $seq->trunc(@$keep)->revcom->seq,
		    ),"\n";
	    }
	}	
    }
}
