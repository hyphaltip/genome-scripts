#!/usr/bin/perl -w
use strict;

my $infile = shift || "shared_target_cleaned";
open(my $fh => $infile) || die $!;
while(<$fh>) {
    chomp;
    my ($dest,$target) = split(/,/,$_,2);
    next unless $target =~ /stajich/;
#    warn("$target -> $dest\n");
    if( -l $target ) {
	my $l = readlink($target);
	warn("link is $l, dest is $dest\n");
	my $oldest = $l;
	if( $l =~ s!/home_stajichlab/jstajich/bigdata!/bigdata/jstajich!) {
	    warn("-> updating $target to point to $l, instead of $oldest\n");
	    unlink($target);
	    symlink($l,$target);
	} elsif( $l =~ s!/home_stajichlab/jstajich!/rhome/jstajich!) {
	    warn("-> updating $target to point to $l, instead of $oldest\n");
	    unlink($target);
	    symlink($l,$target);
	} elsif( $l =~ s!/srv/zpools/tern.ib_bigdata/home/stajichlab/shared!/shared/stajichlab! ) { 
	    warn("-> updating $target to point to $l, instead of $oldest\n");
	    unlink($target);
	    symlink($l,$target);	    
	} else {
	    warn("  not updating $l ($target)\n");
	}
   } else {
#	warn("$target does not exist or is not a link\n");
    }
}
 
