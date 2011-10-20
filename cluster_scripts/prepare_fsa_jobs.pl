#!/usr/bin/perl -w
use strict;
use File::Spec;
use Env qw(HOME);
my $dir = shift || die $!;
$dir = File::Spec->rel2abs($dir);
my $odir = "$HOME/jobs";
my $prefix = 'fsa';

opendir(DIR, $dir) || die $!;

for my $n ( readdir(DIR) ) { 
    next unless $n =~ /^(\d+)$/;
    my $ct = $1;
    open(my $jobfh => ">$odir/$prefix.$ct.sh") || die $!;
    print $jobfh "#PBS -N $prefix.$ct\n","#PBS -l nodes=1:ppn=1\n";
    print $jobfh "cd $dir/$n\n";
    print $jobfh "fsa --mercator cons --exonerate seqs.fasta > seqs.fsa.mfa\n";
    print $jobfh "ln -s seqs.fsa.mfa mavid.mfa\n";
    print $jobfh "fa2phylip < seqs.fsa.mfa > seqs.fsa.phy\n";
    print $jobfh "ln -s seqs.fsa.phy mavid.phy\n";
    close($jobfh);
}
