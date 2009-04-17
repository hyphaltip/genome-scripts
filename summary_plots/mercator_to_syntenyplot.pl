#!/usr/bin/perl -w
use strict;

my $mapfile = shift || die "need a mercator file";

my $ref_genome_index = 0;
my $lb_genome_index = 4;
my $pc_genome_index = 8;

open(my $fh => $mapfile ) || die $!;

open(my $of1 => ">plot/cc-lb_synteny.tab") || die $!;
open(my $of2 => ">plot/cc-pc_synteny.tab") || die $!;
print $of1 '#', join("\t", qw(chrom src_start src_stop strand lb lb_start lb_stop lb_strand)),"\n";

print $of2 '#', join("\t", qw(chrom src_start src_stop strand pc pc_start pc_stop pc_strand)),"\n";

while(<$fh>) {
    my ($aln_id, @line) = split;
    my ($chrom,$start,$end,$strand) = map { $line[$ref_genome_index + $_] } 0..3; 
    my ($lb_chrom,$lb_start,$lb_end,$lb_strand) = map { $line[$lb_genome_index + $_] } 0..3; 
    my ($pc_chrom,$pc_start,$pc_end,$pc_strand) = map { $line[$pc_genome_index + $_] } 0..3;
    next if $chrom eq 'NA';
    
    if( $lb_chrom ne 'NA' ) {
	print $of1 join("\t", $chrom,$start,$end,$strand, # Cc
			$lb_chrom,$lb_start,$lb_end,$lb_strand, # Lb
			),"\n";
    }
    if( $pc_chrom ne 'NA' ) {
	print $of2 join("\t", $chrom,$start,$end,$strand, # Ci
			$pc_chrom,$pc_start,$pc_end,$pc_strand, # Ur
			),"\n";
    }
}
