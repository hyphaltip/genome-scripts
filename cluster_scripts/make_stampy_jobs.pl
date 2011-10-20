#!/usr/bin/perl -w
use strict;
use Env qw(HOME);
use File::Spec;
use Getopt::Long;

my $jobdir = File::Spec->catdir($HOME,'jobs');
my $reference_genome;
my ($outdir,$genome_pref);
my $force = 0;
$outdir = ".";
my $solexaqual = 0;
my $tmpdir = "$HOME/bigdata/tmp/stampy";

#my $otheropts = '--substitutionrate=0.05';

GetOptions(
    'r|ref:s'    => \$reference_genome,
    'o|outdir:s' => \$outdir,
    'p|prefix:s' => \$genome_pref,
    'sq|solexaqual!' => \$solexaqual,
    't|tmp:s'     => \$tmpdir,
    'j|jobdir:s'  => \$jobdir,
    'f|force!'  => \$force);

$outdir = File::Spec->rel2abs($outdir);
if( ! defined $reference_genome )  {
    die "must provide reference genome with -r or --ref\n";
} 
$reference_genome = File::Spec->rel2abs($reference_genome);
my @fqfiles = @ARGV;
my $opts = $solexaqual ? "--solexa " : "";
for my $fqfile ( @fqfiles ) {  
    $fqfile = File::Spec->rel2abs($fqfile);
    my (undef,$bdir,$fname) = File::Spec->splitpath($fqfile);
    next unless ( $fname =~ /(\S+)\.fq$/ || $fname =~ /(\S+)\.fq\.gz$/);
    my $base = $1;
    $base =~ s/\.trim//;
#    $bdir =~ s/\/$//;
    open(my $jobfh => ">$jobdir/$base.$genome_pref.stampy.sh") || die $!;
    print $jobfh "#PBS -N $base\n";
    print $jobfh "cd $bdir\n";
    if( $force ) {
	unlink("$outdir/$base.$genome_pref.sam");
    }    
    print $jobfh "if [ ! -f $fqfile.recaldata ]; then\n stampy.py -g $reference_genome -h $reference_genome -R $fqfile \nfi\n";    
    print $jobfh "if [ ! -f $outdir/$base.$genome_pref.sam ]; then\n stampy.py $opts -g $reference_genome -h $reference_genome --readgroup=ID:$base,SM:$base --baq --bwaoptions=\"-q10 $reference_genome\" --bwatmpdir=$tmpdir -M $fqfile -f sam -o $outdir/$base.$genome_pref.sam\nfi\n";
}
