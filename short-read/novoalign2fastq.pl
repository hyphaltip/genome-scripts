#!/usr/bin/perl -w
use strict;
use Getopt::Long;

while(<>) {
    next if /^\#/;
    chomp;
    my ($id,$status,$seq,$qual,$matches) = split(/\t/,$_);
    next if ($matches eq 'QC' || length($seq) == 0);
    $id =~ s/>//;
    my $len = length($seq);
    
    print join("\n",$id,$seq,"+",substr($qual,0,$len)),"\n";
}
