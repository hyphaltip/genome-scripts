#!/usr/bin/perl -w
use strict;
use Env qw(HOME);
use File::Spec;
use Getopt::Long;

my $jobdir = File::Spec->catdir($HOME,'jobs');
my $reference_genome;
my $gatk_jarfile = '/srv/zpools/tern.ib_bigdata/home/stajichlab/shared/pkg/GATK/GenomeAnalysisTK.jar';
my $javamem ="-Xmx1g";
my $tmpdir  = " -Djava.io.tmpdir=~/bigdata/tmp";
my $force = 0;
GetOptions(
    'r|ref:s'    => \$reference_genome,
    'jar|gatk:s' => \$gatk_jarfile,
    'f|force!'  => \$force);

if( ! defined $reference_genome )  {
    die "must provide reference genome with -r or --ref\n";
} 
$reference_genome = File::Spec->rel2abs($reference_genome);
my @bamfiles = @ARGV;
for my $bam ( @bamfiles ) {  
    $bam = File::Spec->rel2abs($bam);
    my (undef,$bdir,$fname) = File::Spec->splitpath($bam);
    next unless ( $fname =~ /(\S+)\.bam$/);
    my $base = $1;
#    $bdir =~ s/\/$//;
    open(my $jobfh => ">$jobdir/$base.gatk_realign.sh") || die $!;
    print $jobfh "#PBS -N $base\n";
    print $jobfh "cd $bdir\n";
    if( $force ) {
	unlink("$bdir/$base.intervals");
    }
    print $jobfh "if [ ! -f $base.intervals ]; then\n  java $tmpdir $javamem -jar $gatk_jarfile -T RealignerTargetCreator -R $reference_genome -o $base.intervals -I $fname; fi\n";

    if( $force ) {
	unlink("$bdir/$base.realign.bam");
    }    
    print $jobfh "if [ ! -f $base.realign.bam ]; then\n  java $tmpdir $javamem -jar $gatk_jarfile -T IndelRealigner -R $reference_genome -targetIntervals $base.intervals -I $fname -o $base.realign.bam; fi\n";
    
    print $jobfh "if [ ! -f $base.snps.raw.vcf ]; then\n  java $tmpdir $javamem -jar $gatk_jarfile -T UnifiedGenotyper -R $reference_genome -I $base.realign.bam -o $base.snps.raw.vcf -glm BOTH -rf BadCigar ; fi\n";
}
