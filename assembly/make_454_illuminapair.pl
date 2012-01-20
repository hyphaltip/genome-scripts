#!/usr/bin/perl -w
use strict;

while(<>) {
    my $id = $_;
    my $seq = <>;
    my $dsc = <>;
    my $qual = <>;
    if( $id =~ s/\.[fF]/\/1/ ) {
    } elsif( $id =~ s/\.[rR]/\/2/ ) {
    } else {
	warn("no match\n");
    }
    print $id,$seq,$desc,$qual;
}
    
