#!/usr/bin/perl
use warnings;
use strict;
use File::Spec;

my $jobdir = 'jobs';
my $genome = '/home_stajichlab/jstajich/bigdata/Bd/genomes/JEL423_17Jan07.cs.fa';
my $indexes = '1,2,3,4,5,6,7,8,9,10';
my $threads = '4';
my $space = '-A 1';
my $match_opts = '-l -T /dev/shm/ -M 3 -w 0';
my $postproc_opts = "-O 1 -S 0 -a 2 -Y 2 -P 2";
my $odir = "bfast_aln";

for my $full( @ARGV ) {
 my (undef,$dir,$file) = File::Spec->splitpath($full);
 warn("file is $file\n");
 next unless( $file =~ /(\S+)\.csfastq/);
 my $stem = $1;
 open(my $fh => ">$jobdir/$stem.bfast.JEL423.sh") || die $!;
 print $fh "#PBS -N $stem.bfast -q highmem -l nodes=1:ppn=$threads,walltime=99:99:99,mem=6gb -j oe\n";
 print $fh "cd /home_stajichlab/jstajich/bigdata/Bd/April_2012/SOLiD\n";
 print $fh "module load stajichlab\n","module load bfast\n";
 print $fh "if [ ! -f $odir/$stem.JEL423.bmf ]; then \n";
 print $fh "bfast match -f $genome -n $threads $space $match_opts -r $full > $odir/$stem.JEL423.bmf\n";
 print $fh "fi\n";
 print $fh "if [ ! -f $odir/$stem.JEL423.baf ]; then\n";
 print $fh "bfast localalign -f $genome -n $threads $space -m $odir/$stem.JEL423.bmf > $odir/$stem.JEL423.baf\n";
 print $fh "fi\n";
 print $fh "if [ ! -f $odir/$stem.JEL423.sam ]; then\n";
 print $fh "bfast postprocess -f $genome -n $threads $space $postproc_opts -i $odir/$stem.JEL423.baf > $odir/$stem.JEL423.sam\n";
 print $fh "fi\n";
}
