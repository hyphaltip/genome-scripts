#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $groups;
GetOptions(
	   'g|lookup:s' => \$groups);

my %groups;
if( $groups ) {
    open(my $fh => $groups) || die "cannot open $groups: $!";
    while(<$fh>) {
	next if /^#/;
	my ($grpname,@members) = split;
	for my $m ( @members ) { $groups{$m} = $grpname }
    }
}
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
    my $str;
    if( keys %count > 1 ) {
	my %r;
	for my $k ( keys %count ) {
	    if(exists $groups{$k} ) {
		$r{$groups{$k}}++;
	    }
	}	
	$str =  join("-", sort keys %r);	
    } else {
	$str =  join("-", sort keys %count);
    }
    push @{$patterns{$str}}, $group;
}

for my $p ( sort { scalar @{$patterns{$a}} <=>
		       scalar @{$patterns{$b}} }  keys %patterns ) {
    print join("\t",$p,scalar @{$patterns{$p}}),"\n";
    open(my $ofh => ">$file.patterns/$p.dat") || die $!;
    print $ofh join("\n", @{$patterns{$p}}),"\n";
}

