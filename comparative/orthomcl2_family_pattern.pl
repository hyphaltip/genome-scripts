#!/usr/bin/perl -w
use strict;

my $file = shift;
mkdir("$file.patterns");
my %sp;
my %patterns;
open(my $fh => $file) || die $!;
while(<$fh>) {
    my ($group,@orthologs) = split;
    $group =~ s/://;
    my %count;
    my %domains;	
    for my $gene ( @orthologs ) {
	my ($sp)= split(/\|/,$gene);
	$count{$sp}++;
	$sp{$sp} = 1; # collecting all the species names
    }
    my $str = join("-", sort keys %count);
    push @{$patterns{$str}}, $group;
}

for my $p ( sort { scalar @{$patterns{$a}} <=> scalar @{$patterns{$b}} }  keys %patterns ) {
    print join("\t",$p,scalar @{$patterns{$p}}),"\n";
    open(my $ofh => ">$file.patterns/$p.dat") || die $!;
    print $ofh join("\n", @{$patterns{$p}}),"\n";
}

