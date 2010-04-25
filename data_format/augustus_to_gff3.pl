#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $prefix = 'augustus';

GetOptions(
	   'p|prefix=s' => \$prefix,
	   );
print "#gff-version 3\n";
while(<>) {
    next if /^\#/;
    next if /^\s+$/;
    next if /^(Make|Fehler|Error)/;
    chomp;
    my @line = split(/\t/,$_);
    if( ! defined $line[2] ) { 
	warn("'$_'\n");
	next;
    }
    if( $line[2] eq 'CDS' || $line[2] eq 'transcript' || $line[2] eq 'gene' ) {
	my $lastcol = pop @line;
	$line[2] = 'mRNA' if $line[2] eq 'transcript';
        $lastcol =~ s/(ID|Parent)=/$1=$prefix-/g;
	if( $line[2] eq 'CDS') {
		$lastcol =~ s/ID=([^;]+);//;
	} else {
		if( $lastcol =~ /ID=([^;]+);?/ ) {
		  $lastcol .= ";Name=$1";
		}
	}
	print join("\t", @line, $lastcol),"\n";
    }
}
