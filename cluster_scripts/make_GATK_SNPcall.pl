#!/usr/bin/perl -w
use strict;
use Env qw(HOME);
use File::Spec;
use Getopt::Long;

my $jobdir = File::Spec->catdir($HOME,'jobs');
my $reference_genome;
my $gatk_jarfile = '/srv/zpools/tern.ib_bigdata/home/stajichlab/shared/pkg/GATK/GenomeAnalysisTK.jar';
my $javamem ="-Xmx16";
my $tmpdir  = " -Djava.io.tmpdir=~/bigdata/tmp";
my $force = 0;
my $outname = '../vcf/allSNPS.raw.vcf';
my $outmetrics = '../vcf/allSNPs.info';
GetOptions(
    'r|ref:s'    => \$reference_genome,
    'jar|gatk:s' => \$gatk_jarfile,
    'o|out:s'    => \$outname,
    'f|force!'  => \$force);

if( ! defined $reference_genome )  {
    die "must provide reference genome with -r or --ref\n";
} 
print "#PBS -N GATKsnps\n";
$reference_genome = File::Spec->rel2abs($reference_genome);
my @bamfiles = @ARGV;
my ($first) = @bamfiles;
my (undef,$bdir,$fname) = File::Spec->splitpath(File::Spec->rel2abs($first));
print "cd $bdir\n";    
my $bamlst = join(" ", map { sprintf("-I %s",$_) } @bamfiles);
print "java $tmpdir $javamem -jar $gatk_jarfile -T UnifiedGenotyper -R $reference_genome $bamlst -o $outname -glm BOTH --read_filter BadCigar -nt 48 --metrics_file $outmetrics\n";
