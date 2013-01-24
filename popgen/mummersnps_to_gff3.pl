#!/usr/bin/perl
use warnings;
use strict;

my $source = 'MUMMER';
my $snptype   = 'SNP';
my $indeltype = 'INDEL';

my %count = ($snptype => 1,
	     $indeltype => 1);
<>;
my $hdr = <>;
if( $hdr !~ /^NUCMER/) {
    warn "expected NUMBER as $hdr\n";
}
while(<>) {
    next if /^(NUCMER|\s+|\[)/;
    my @row = split;
    my $type = $snptype;
    if( $row[1] eq '.' || $row[2] eq '.') {
	$type = $indeltype;
    }
    print join("\t", $row[10],$source, $type, 
	       $row[0], $row[0]+1,
	       '.',
	       '.',
	       '.',
	       join(";",
		    sprintf("ID=%s%07d",$type,$count{$type}++),
		    sprintf("Target=%s %d %d %s",$row[11], $row[3],
			    $row[3] + 1,
			    $row[9] > 0 ? '+' : '-'))),"\n";
}
