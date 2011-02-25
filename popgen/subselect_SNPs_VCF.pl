#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use List::Util qw(sum shuffle);
use Getopt::Long;
my $target_freq_min = 0.40;
my $target_freq_max = 0.60;
my $flank = 50;

my (%locs,%dat);
my $total_SNPs;
my @strains;

GetOptions(
	   'min|minfreq:s' => \$target_freq_min,
	   'max|maxfreq:s' => \$target_freq_max,
	   'f|flank:i'     => \$flank,
	   
	   );
	   
my ($ref,$vcf) = @ARGV;

my $db = Bio::DB::Fasta->new($ref);

open(my $fh => $vcf) || die $!;

while(<$fh>) {
    if( /^\#CHROM/) {
	my @header = split;
	splice(@header,0,9);
	@strains = @header;
    }
    next if /^\#/;
    my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@data) = split;
    next if length($alt) > 1;
    next if $info =~ /^INDEL/;
    my %inf = map { split(/=/,$_,2) } split(/;/,$info);
    next if $inf{AF1} < $target_freq_min || $inf{AF1} > $target_freq_max;
    push @{$dat{$chrom}}, [$pos,$id,$ref,$alt,\%inf,\@data];
    $locs{$chrom}->{$pos}++;
    $total_SNPs++;
}
warn "$total_SNPs SNPs count\n";

print join("\t", qw(CHROM POS REF_ALLELE ALT_ALLELE ML_ALLELE_FREQ 
		    ALT_ALLELE_COUNT LEFT_FLANK REF_BASE RIGHT_FLANK)), "\n";

for my $chrom ( sort keys %dat ) {
    my $ct = 0;
    for my $p ( shuffle @{$dat{$chrom}} ) {
	my $pos = $p->[0];
	my ($flank_left_5, $flank_left_3) = ($pos-$flank,$pos-1);
	my ($flank_right_5, $flank_right_3) = ($pos+1,$pos+$flank);
	
	my $skip = 0;
	for( my $i = $flank_left_5; $i <= $flank_left_3; $i++ ) {
	    if( exists $locs{$chrom}->{$i} ) {
		$skip = 1;
		last;
	    }
	}
	next if $skip;

	for( my $i = $flank_right_5; $i <= $flank_right_3; $i++ ) {
	    if( exists $locs{$chrom}->{$i} ) {
		$skip = 1;
		last;
	    }
	}
	next if $skip;
	next if $flank_left_5 < 1 || $flank_right_3 > $db->length($chrom);
	my $vals = $p->[5];
	my @allelecodes = map { (split(':',$_))[0] } @$vals;
	my $sum = sum (map { split(/\//,$_) } @allelecodes);
	print join("\t", $chrom, $pos, $p->[2], $p->[3],
		   sprintf("AF1=%.4f",$p->[4]->{AF1}),
		   $sum,
		   $db->seq($chrom,$flank_left_5 => $flank_left_3),
		   lc $db->seq($chrom,$pos => $pos),		   
		   $db->seq($chrom,$flank_right_5 => $flank_right_3)),"\n";
	last if $ct++ > 10;
    }
}
