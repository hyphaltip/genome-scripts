#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;
my $jobdir = 'jobs';
my $ref = 'JEL423_17Jan07.fa';
my $fname = "mpileup";
GetOptions('j|jobs:s' => \$jobdir,	   
	   'r|ref:s'  => \$ref,
	   'f|fname:s' => \$fname,
    );

mkdir($jobdir);
$ref = File::Spec->rel2abs($ref);

for my $file ( @ARGV ) {
    $file = File::Spec->rel2abs($file);
}

#    my (undef,undef,$fname) = File::Spec->splitpath($f);
#    $fname =~ s/\.bam$//;
#    print "$f\n";

open(my $jobfh => ">$jobdir/$fname.sh") || die $!;
print $jobfh "#PBS -N $fname -q highmem -j oe -l nodes=1:ppn=1,mem=2gb,walltime=48:00:00\n";
print $jobfh join("\n", map { sprintf("module load %s",$_) }
		  qw(stajichlab samtools varscan java) ),"\n";
print $jobfh "samtools mpileup -f $ref @ARGV > $fname.mpileup\n";
print $jobfh "#java -jar \$VARSCANJAR mpileup2snp --output-vcf > $fname.varscan.vcf\n";

#warn("here\n");
