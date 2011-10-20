#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use Env qw(HOME);
my $jobdir = "$HOME/jobs";
my $odir   = "$HOME/output";
my $exe = 'hmmscan-3.0';
my $dir;
my $db = '/srv/projects/db/pfam/2009-10-Pfam24.0/Pfam-A.hmm';
my $cpus = 4;
GetOptions(
	   'db:s'    => \$db,
	   'd|dir:s' => \$dir);

$dir = File::Spec->rel2abs($dir);
opendir(DIR, $dir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.(fa|fasta|pep|seq)$/;
    my $stem = $1;
    open(JOB, ">$jobdir/$stem.pfamscan.sh") || die $!;
    print JOB "#PBS -l nodes=1:ppn=$cpus\n";
    print JOB "$exe --cpu $cpus --tblout $odir/$stem.tbl --domtblout $odir/$stem.domtbl.tab $db $dir/$file > $odir/$stem.Pfam_24.hmmscan\n";
    close(JOB);
}
