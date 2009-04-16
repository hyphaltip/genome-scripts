#!/usr/bin/perl -w
# $Id: rptmask_to_summary.pl 250 2008-12-12 22:20:20Z stajich $
use strict;


print join("\t",qw(#chrom
		   start stop class)),"\n";
while(<>) {
    next if( /^(\s+SW|score)/ || /^\s+$/);
    my ($swscore,$perdiv,$perdel,$perins,
	$qid,$qstart,$qend,$left, $strand,
	$hid, $class, $hstart,$hend, $hleft, $uid) = split;
    if( $class eq 'Low_complexity' || $class eq 'Simple_repeat') {
	print join("\t", $qid, $qstart,$qend, $class ),"\n";
    } else {
 	print join("\t", $qid, $qstart,$qend, "$class:$hid"),"\n";
    }
}
