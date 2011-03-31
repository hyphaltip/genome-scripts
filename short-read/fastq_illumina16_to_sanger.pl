#!/usr/bin/perl -w
use strict;

my $state = 0;
while(<>) {
    if( /^@/ ) {
	print;
	$state=1;
    }elsif( $state == 1 ) {
	print ;
	$state++; 
    } elsif( $state == 2 ) {
	print "+\n";
	$state++;
    } elsif( $state == 3) {
	chomp;
	my $qual = $_;
	my $len = length($qual);
	for( my $i = 0; $i < $len; $i++ ) {
	    my $char = substr($qual,$i,1);
	    my $old_val = ord($char);
	    my $new_val = $old_val - 64 + 33;	    
#	    warn("$char => $old_val -> $new_val\n");
	    substr($qual,$i,1,chr($new_val));
	}
	print $qual,"\n";
#	die;
	$state++;
    } else {
	die("unknown state $state with $_");
    }
}
