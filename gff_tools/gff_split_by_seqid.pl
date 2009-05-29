#!/usr/bin/perl -w
use strict;

my %seqids;

while(<>) {
    my ($seqid) = split(/\t/,$_);
    push @{$seqids{$seqid}}, $_;
}
while( my ($seqid,$dat) = each %seqids ) {
    open(my $fh => ">$seqid.gff") || die $!;
    for my $r ( @{$dat} ) {
	print $fh $r;
    }
}
    
