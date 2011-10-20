#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use Cwd;
my $db;

my $cwd = getcwd;
GetOptions('d|db=s' => \$db);

my $fullname = File::Spec->rel2abs($db);
my (undef,undef,$dbname) = File::Spec->splitpath($fullname);
# remaining file names are the query files
for my $q ( @ARGV) {
    my $stem = $q;
    $stem =~ s/\.(fa|fasta|fastq|fast|dna|seq|txt|fq)$//;
    
    my $base = "$stem.$dbname";
    open(JOB, ">$base.sh")|| die $!;
    print JOB "#!/bin/bash\n",    
    "#PBS -l nodes=1:ppn=1\n",
    "#PBS -N $stem\n",
    "cd $cwd\n",
    "bwa aln -f $base.sai $fullname $q\n",
    "bwa samse $fullname $base.sai $q | gzip -c > $base.sam.gz\n",
    "samtools import ITAG1_genomic.fai -o $base.bam $base.sam.gz\n",
    "samtools sort $base.bam\n",
    "samtools index $base.bam\n";
    close(JOB);
}
