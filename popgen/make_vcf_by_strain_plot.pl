#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
my %dat;
my @strains;
while(<>) {
    if( /^\#CHROM/) {
	my @header = split;
	splice(@header,0,9);
	@strains = @header;
    }
    next if /^\#/;
    my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@data) = split;
    push @{$dat{$chrom}}, [$pos,$id,$ref,$alt,$info,\@data];
}
for my $strain ( @strains ) {
    my $i =0;
    open(my $fh => ">$strain.var.dat") || die $!;
    for my $chrom ( sort keys %dat ) {
	for my $p ( @{$dat{$chrom}} ) {
	    my $vals = $p->[5];
	    my ($allelecode) = split(':',$vals->[$i]);
	    my $sum = sum split(/\//,$allelecode);
	    print $fh join("\t", $chrom, $p->[0],$sum / 2.0);
	}
    }
    $i++;
}
