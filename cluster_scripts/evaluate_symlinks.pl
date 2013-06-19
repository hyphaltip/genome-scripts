#!/usr/bin/perl -w
use strict;

while(<>) {
    my ($link) = split;
    my $readlink = readlink($link);
    if( $readlink =~ /zfs/ ||
	$readlink =~ /home_stajichlab/ ) {
	print "$link,$readlink\n";
    }
}
