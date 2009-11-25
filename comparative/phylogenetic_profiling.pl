#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SearchIO;
use File::Spec;
use Bio::DB::Fasta;

my $tab_dir; # directory of tab delimited
my $ext = 'tab';
my $db;
my $paligned = 0.20;
my $cutoff_evalue = '0.001';
my $odir = '.';
GetOptions(
	   't|d|tab:s' => \$tab_dir,
	   'o|output:s' => \$odir,
	   'e|evalue:s' => \$cutoff_evalue,
	   );

opendir(DIR, $tab_dir) || die $!;

my %results;
my %hsp_list;
my %qsp_list;
for my $f ( readdir(DIR) ) {
    next unless $f =~ /\.$ext$/;
    open(my $fh => File::Spec->catfile($tab_dir,$f)) || die "$f: $!";
    while(<$fh>) {
	next if /^\#/;
	next if /^\s+$/;
	my ($qname,$hname, $percent_id, $hsp_len, $mismatches,$gapsm,
	    $qstart,$qend,$hstart,$hend,$evalue,$bits) = split;
	next if $evalue > $cutoff_evalue;
	my ($qsp,$qseq) = split('_',$qname,2);
	my ($hsp,$hseq) = split('_',$hname,2);
	$hsp_list{$hsp}++;
	$qsp_list{$qsp}++;
	push @{$results{$qsp}->{$qname}->{$hsp}->{$hname}}, [$percent_id,$hsp_len,$qstart,$qend,$hstart,$hend,$evalue,$bits];
    }
}
my @qsp = sort keys %qsp_list;
my @hsp = sort keys %hsp_list;
for my $qsp ( @qsp ) {
    open(my $ofh => ">$odir/$qsp.phyloprofile.tab") || die $!;
    print $ofh join("\t", 'GENE', @hsp), "\n";
    for my $gene ( sort keys %{$results{$qsp} || {}} ) {
	print $ofh join("\t", $gene, map { sprintf("%d",(scalar keys %{$results{$qsp}->{$gene}->{$_} || {}})||0) } @hsp ),"\n";
    }
}
