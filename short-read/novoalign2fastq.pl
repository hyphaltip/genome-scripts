#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $code = '200M6AAXX_4_R';

GetOptions('c|code:s' => \$code);
while(<>) {
    next if /^\#/;
    chomp;
    my ($id,$status,$seq,$qual,$matches) = split(/\t/,$_);
    next if ($matches eq 'QC' || length($seq) == 0);
    $id =~ s/>//;
    my $len = length($seq);
    
    print join("\n",'@'.$code.$id,$seq,"+",substr($qual,0,$len)),"\n";
}
