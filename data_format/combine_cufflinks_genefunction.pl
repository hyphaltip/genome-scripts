#!/usr/bin/perl -w
use strict;
use warnings;

my ($fasta,$cuff) = @ARGV;


my %lookup;
open(my $fh => $fasta) || die $!;
while(<$fh>) {
    if( />(\S+).+Name:\"([^"]+)\"/ ) {
	my ($id,$name) = ($1,$2);
	$lookup{$id} = $name;
    }
}
open($fh=> $cuff) || die $!;

my $header =<$fh>;
chomp($header);
print $header,"\t","Function","\n";
while(<$fh>) {
    my @row = split;
    push @row, $lookup{$row[0]};
    print join("\t", @row), "\n";
}

