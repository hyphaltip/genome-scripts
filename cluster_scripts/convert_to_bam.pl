#!/usr/bin/perl -w
use strict;
use warnings;
use Env qw(HOME);
use File::Spec;
use Getopt::Long;

my $jobdir = File::Spec->catdir($HOME,'jobs');
my $force = 0;
my $outdir = ".";

GetOptions(
    'o|outdir:s' => \$outdir,
    'j|jobdir:s'  => \$jobdir,
    'f|force!'  => \$force);

$outdir = File::Spec->rel2abs($outdir);
my @samfiles = @ARGV;
for my $samfile ( @samfiles ) {  
    $samfile = File::Spec->rel2abs($samfile);
    my (undef,$bdir,$fname) = File::Spec->splitpath($samfile);
     warn("fname is $fname\n");
    next unless ( $fname =~ /(\S+)\.sam$/);
    my $base = $1;
    open(my $jobfh => ">$jobdir/$base.makebam.sh") || die $!;
    print $jobfh "#PBS -N $base.makebam\n";
    print $jobfh "cd $bdir\n";
    unless( $force ) {
	print $jobfh "if [ ! -f $base.bam ]; then\n";    
    } 
#    print $jobfh "  java -jar $HOME/bigdata-shared/pkg/Picard/picard/SamFormatConverter.jar INPUT=$fname OUTPUT=$base.unsrt.bam\n";
 #   print $jobfh "  java -jar $HOME/bigdata-shared/pkg/Picard/picard/SortSam.jar INPUT=$base.unsrt.bam OUTPUT=$base.bam SORT_ORDER=coordinate\n";
    print $jobfh "  samtools view -b -S $fname > /scratch/$base.unsrt.bam\n";
#    print $jobfh "  java -jar $HOME/bigdata-shared/pkg/Picard/picard/SortSam.jar INPUT=/dev/shm/$base.unsrt.bam OUTPUT=$base.bam SORT_ORDER=coordinate\n";
    print $jobfh "  samtools sort /scratch/$base.unsrt.bam /scratch/$base\n";
    print $jobfh "  mv /scratch/$base.bam . \n";
    print $jobfh "  /bin/rm -f /scratch/$base.unsrt.bam\n";
    unless( $force ) {
	print $jobfh "fi\n";
    }
}
