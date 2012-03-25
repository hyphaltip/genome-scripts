#!/usr/bin/perl -w
use strict;

open(my $skipfh => ">skip") || die $!;

while(my $id = <>) {
    my $seq = <>;
    my $desc =<>;
    my $qual = <>;
    
    if (($id =~ /\S+\s+\d+\:([YN]):\d+/) && $1 eq 'N') {
	print $skipfh $id,$seq,$desc,$qual;
	next;
    }
    
    print $id,$seq,$desc,$qual;
}
