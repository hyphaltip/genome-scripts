#!/usr/bin/perl -w
use strict;
use Bio::Root::IO;
use Getopt::Long;

my ($run_wise,$run_fasta, $run_exonerate) = (0,1,1);

GetOptions(
	   'w|wise!' => \$run_wise,
	   'f|fasta|fastx!' => \$run_fasta,
	   'e|exonerate!' => \$run_exonerate);
	   
my $io = Bio::Root::IO->new;
my $pseudowise = $io->exists_exe('pseudowise');
my $tfastx = $io->exists_exe('tfastx36_t');
my $exonerate = $io->exists_exe('exonerate');

my $dir = shift || die $!;

opendir(my $dirh => $dir) || die $!;
for my $file ( readdir($dirh) ) {
    next unless ($file =~ /(\S+)\.genomic.fa/);
    my $stem = $1;
    warn("$stem\n");
    
    `$pseudowise -pretty -quiet $dir/$1.protein.fa $dir/$1.cds.fa $dir/$file >$dir/$1.pseudowise` if $run_wise && $pseudowise;
    `$tfastx -H -E 1e-2 $dir/$1.protein.fa $dir/$file > $dir/$1.TFASTX` if $run_fasta && $tfastx;
    `$exonerate -f -10 --refine full --showtargetgff -m protein2genome -q $dir/$1.protein.fa -t $dir/$file >$dir/$1.EXONERATE` if $run_exonerate && $exonerate;
}

