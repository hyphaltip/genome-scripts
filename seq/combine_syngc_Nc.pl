#!/usr/bin/perl -w
use strict;

my ($chips,$gc) = @ARGV;

open(my $fh => $chips) || die $!;
my %cds;
while(<$fh>) {
    my ($gene, undef,undef,$Nc) = split;
    $cds{$gene}->{Nc} = $Nc;
}
open($fh => $gc) || die $!;
while(<$fh>) {
    my ($gene,$gc) = split;
    $gene =~ s/^([^:]+)://;
    if( ! exists $cds{$gene} ) {
	warn("no $gene in chips\n");
	next;
    }
    $cds{$gene}->{GC} = $gc;
}
print join("\t", qw(GENE CHIPS GC)),"\n";
for my $gene ( sort keys %cds ) {
  next unless exists $cds{$gene}->{Nc} && $cds{$gene}->{GC};
    print join("\t", $gene, $cds{$gene}->{Nc}, $cds{$gene}->{GC}),"\n";
}
