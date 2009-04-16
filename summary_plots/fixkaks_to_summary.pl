#!/usr/bin/perl -w
use strict;

# Hardcode for Tom S's file from the cocci project
my $genefile = 'plot/ci_gene_summary.tab';
open(my $fh => $genefile) || die $!;
my %genes;
my $header = <$fh>;
$header =~ s/^\#//;
my $i =0;
my %header = map { $_ => $i++ } split(/\t/,$header);
my %gene;
while(<$fh>) {
    my @row = split;
    $gene{$row[ $header{'gene'} ]} = [$row[ $header{'chrom'} ],
				      $row[ $header{'chrom_start'} ]];
}

my @final;
while(<>) {
    chomp;
    my ($immi,$posa,$ka,$ks,$omega) = split(/\t/,$_);
    if( $immi =~ /^\d+\.m\d+/ ) {
	($immi,$posa) = ($posa,$immi);
    }
    if( $immi !~ s/cdna_// ) {
	warn("C.immitis gene name ($immi) is not what we expected\n");
    }
    $immi .= ".2"; # force to verion 2 
    if( ! exists $gene{$immi} ) {
	warn("cannot find gene $immi from summary file\n");
	next;
    }
    push @final, [ $gene{$immi}->[0], $immi, $posa, $gene{$immi}->[1],
		   $ka,$ks,$omega];
}

print join("\t",'#chrom',qw(cimmitis cposadasii start ka ks kaks)),"\n";

for my $f ( sort { $a->[0] cmp $b->[0] ||
		       $a->[3] <=> $b->[3] }
	    @final ) {
    print join("\t", @$f), "\n";

}
