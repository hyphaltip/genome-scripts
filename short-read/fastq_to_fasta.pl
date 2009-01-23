#!/usr/bin/perl -w
use strict;
my $ready = 0;
while(<>) {    
    if($ready == 0 ) {
	if( /^@(\S+)/) {
	    print ">$1\n";
	} else {
	    die("off track $_");
	}
    } elsif( $ready == 1 ) {
	print;
    } 
    $ready++;
    $ready = 0 if $ready == 4;
}

