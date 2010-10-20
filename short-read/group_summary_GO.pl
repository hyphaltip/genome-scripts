#!/usr/bin/perl -w
use strict;
my $cutoff_evalue = 1e-5;

my $gene2pfam = '/project/tomato/ITAG_1/10-March-2010/ITAG1_proteins.Pfam_24.tbl';
my $pfam2go = '/project/db/GO/pfam2go';

open(my $gfh => $gene2pfam) || die $!;
my %gene2pfam;
while(<$gfh> ){  # hmmer3 .tbl table
    next if /^\#/;
    my ($domain,$acc,$gene,$geneacc,$evalue) = split(/\s+/,$_);
    next if $evalue > $cutoff_evalue;
    $gene2pfam{$gene}->{$domain}++;    
}
close($gfh);
open($gfh => $pfam2go) || die $!;
my %pfam2go;
while(<$gfh>) {
    next if /^\!/;
    chomp;
    my ($pfam,$desc) = split(/\s+>\s+/,$_);
    my ($pfamacc,$pfamname) = split(/\s+/,$pfam);
    my ($goterm, $goid) = split(/\s+;\s+/,$desc);
    
    push @{$pfam2go{$pfamname}}, [$goid,$goterm,$pfamname];
    
}
my $dir = shift || ".";

opendir(DIR, $dir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.dat$/;
    my $stem = $1;
    open(my $fh => "$dir/$file") || die $!;
    my %go_seen;
    my %genepfam;    
    while(<$fh>) {
	chomp;
	next if /^FNAME/;
	my ($gene,@exp) = split(/\t/,$_);
	my $desc = pop @exp;
	
	$gene =~ s/^mRNA://;

	for my $domain ( keys %{$gene2pfam{$gene}||{}} ) {	    
	    for my $go ( @{$pfam2go{$domain} || []} ) {
		my ($goid,$goterm,$pfam) = @$go;
		$go_seen{"$goid\t$goterm"}->{$gene}++;
		$genepfam{$gene}->{$pfam}++;
	    }
	}
    }
    open(my $rpt => ">$stem.GO") || die $!;
    for my $go ( sort { scalar keys %{$go_seen{$b}} <=> scalar keys %{$go_seen{$a}} } keys %go_seen ) {
	my @genes = sort keys %{$go_seen{$go}};
	my %tot;
	for my $k ( map { keys %{$genepfam{$_}} } @genes) { $tot{$k}++ }
	print $rpt join("\t", $go, scalar @genes,
			join(",",map { sprintf("%s=%d",$_,$tot{$_}) } 
			     sort { $tot{$b} <=> $tot{$a} } keys %tot)),"\n";
    }
    close($rpt);
}
