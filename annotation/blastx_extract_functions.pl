#!/usr/bin/perl

# for a metagenome project, extract the BLASTX hits
# matches these against DB of GO/pathways e.g. one produced by interpro2qiime_mds.pl
# prints out a table of the pathways seen and their frequencies
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
my $min_len = 15;
my $ext = '.blast.bl9';
my $dir = shift || die "need the name of a dir with the BLAST results";
my $dbdir = shift || "/srv/projects/db/QIIME/fungal-012013/per-genome/";

opendir(IN, $dir) || die $!;
my %files;
for my $file ( readdir(IN) ) { 
    next unless $file =~ /\Q$ext\E$/;
    my ($pref) = split(/\./,$file);
    push @{$files{$pref}}, File::Spec->catfile($dir,$file);    
}

opendir(DB,$dbdir)|| die $!;
my %gene2function;
for my $f ( readdir(DB) ) {
    next unless $f =~ /\.txt$/;
    open(my $fh => File::Spec->catfile($dbdir,$f)) || die $!;
    while(<$fh>) {
	chomp;
	my ($gene,@functions) = split(/\t/,$_);
	push @{$gene2function{$gene}}, @functions;
    }
}


for my $pref ( keys %files ) {
    my %seenpathways;
    open(my $ofhf => "| gzip -c > $pref.seen.genes.gz") || die $!;
    for my $file ( @{$files{$pref}} ) {
	open(my $fh => $file) || die $!;
	my %hits;
	while(<$fh>) {
	    next if /^\#/;
	    my ($q,$h,$pid,$aln_len,@rest) = split;
	    next unless $aln_len > $min_len;
	    $hits{$q}->{$h} = $rest[-2];
	}
	for my $q ( keys %hits ) {
	    my $seen;
	    print $ofhf join("\t", $q, keys %{$hits{$q}}), "\n";
	    for my $h ( keys %{$hits{$q}} ) {

		if( $gene2function{$h} ) {
		    for my $p ( @{$gene2function{$h}} ) {
			$seenpathways{$p}++;
		    }
		    $seen = 1;
		    last;
		}
	    }
	    if( ! $seen ) {
		$seenpathways{'Unknown_function'}++;
	    }
	}
    }
    open(my $ofh => ">$pref.seen.table") || die $!;
    for my $pth ( sort { $seenpathways{$b} <=> $seenpathways{$a} } 
		       keys %seenpathways ) {
	print $ofh join("\t", $seenpathways{$pth}, $pth), "\n";
    }
}
