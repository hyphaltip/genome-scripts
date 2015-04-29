#!env perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Fasta;
my ($vcf,$db);

GetOptions(
    'i|input|vcf|variant:s' => \$vcf,
    'd|db:s'                => \$db,
    );

my $dbh= Bio::DB::Fasta->new($db);

my $in;
if( $vcf =~ /\.gz/ ) {
    open($in => "zcat $vcf |") || die $!;
} else {
    open($in => $vcf) || die $!;
}

my %counts;
while(<$in>){
    next if /^\#/;
    chomp;
    my ($contig,$pos,$id,$ref,$allele,$qual,$filter,$info) = split(/\t/,$_);
    next if $filter !~ /PASS/;
    $counts{$contig}++;
}

print join("\t", qw(CONTIG TOTAL_SNPS LENGTH SNP_PER_BASE)), "\n";

for my $ctg ( sort { $counts{$b} <=> $counts{$a} } keys %counts ) {
    my $len = $dbh->length($ctg);
    print join("\t", $ctg, $counts{$ctg}, $len, 
	       sprintf("%.2f",1000 * $counts{$ctg}/$len)), "\n";
}
