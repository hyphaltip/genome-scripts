#!/usr/bin/perl -w
use strict;
use Env qw(HOME);
use File::Spec;
use Getopt::Long;

my $jobdir = File::Spec->catdir($HOME,'jobs');
my $reference_genome;
my $gatk_jarfile = '/opt/stajichlab/GATK/latest/GenomeAnalysisTK.jar';
my $javamem ="-Xmx256g";
my $tmpdir  = " -Djava.io.tmpdir=/dev/shm";
my $force = 0;
my $outname = 'allSNPS.raw.vcf';
my $outindel = 'allINDELS.raw.vcf';
my $outmetrics = 'allSNPs.info';
my $outmetricsindel = 'allINDELS.info';
GetOptions(
    'r|ref:s'    => \$reference_genome,
    'jar|gatk:s' => \$gatk_jarfile,
    'o|out:s'    => \$outname,
    'outindel:s' => \$outindel,
    'f|force!'  => \$force);

if( ! defined $reference_genome )  {
    die "must provide reference genome with -r or --ref\n";
} 
open(my $jobfh => ">GATK_run_SNP.sh") || die $!;
print $jobfh "#PBS -N GATKsnps -l nodes=1:ppn=48,walltime=99:99:99,mem=256gb\n","module load stajichlab\n","module load java\n";
$reference_genome = File::Spec->rel2abs($reference_genome);
my @bamfiles = map { File::Spec->rel2abs($_) } @ARGV;
my ($first) = @bamfiles;
my $bamlst = join(" ", map { sprintf("-I %s",$_) } @bamfiles);
print $jobfh "java $tmpdir $javamem -jar $gatk_jarfile -T UnifiedGenotyper -R $reference_genome $bamlst -o $outname -glm SNP --read_filter BadCigar -nt 48 --metrics_file $outmetrics\n";


open($jobfh => ">GATK_run_INDEL.sh") || die $!;

print $jobfh "#PBS -N GATKindel -l nodes=1:ppn=48,walltime=99:99:99,mem=256gb\n","module load stajichlab\n","module load java\n";
print $jobfh "java $tmpdir $javamem -jar $gatk_jarfile -T UnifiedGenotyper -R $reference_genome $bamlst -o $outindel -glm INDEL --read_filter BadCigar -nt 48 --metrics_file $outmetricsindel\n"
