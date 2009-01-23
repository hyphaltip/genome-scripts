#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $extended_eland = 0; # extended eland format

GetOptions(
	   'extended!' => \$extended_eland,
	   );

while(<>) {
    chomp;
    my ($read,$seq,$qual, $no_hits,
	undef, $trimmed_len,
	$strand, $chrom, $start,
	$mm,@mms) = split(/\t/,$_);
    my ($match_type);
    my @hit_dist = (0,0,0);
    if( $mm == 0 ) {
	$hit_dist[0] = $no_hits;
    } elsif( $mm == 1) {
	$hit_dist[1] = $no_hits;
    } elsif( $mm == 2 ) {
	$hit_dist[2] = $no_hits;
    } else {
#	warn("mm = $mm for $_\n");
	next;
    }

    my @rest;
    if( $no_hits == 1 ) {
	$mm = 0;
	$match_type = 'U'.$mm;
	@rest = ($chrom,$start,$strand eq '+' ? 'F' : 'R',
		 "..");
#		 map { my ($refbase,$refloc,$qbase,$qual) = 
#			   ( /([A-Z])\-\>(\d+)([A-Z])(\d+)/ );
#		       sprintf("%d%s",$refloc+1,$refbase);
#		   } @mms);
    } else { 
	$match_type = 'R'.$mm;
    }
    
    print join("\t", ">$read", $seq, $match_type,@hit_dist,@rest), "\n";
	       
}
