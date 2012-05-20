#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use Env qw(HOME);
my $jobdir = "$HOME/jobs";
my $odir   = "$HOME/output";
my $exe = 'hmmscan-3.0';
my $dir = '.';
my $db = '/srv/projects/db/pfam/2011-12-09-Pfam26.0/Pfam-A.hmm';
my $version = "Pfam_26";
my $cpus = 3;
GetOptions(
	   'db:s'    => \$db,
	   'o|odir:s' => \$odir,
	   'version:s' => \$version,
	   'd|dir:s' => \$dir);

$dir = File::Spec->rel2abs($dir);
opendir(DIR, $dir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.(fa|fasta|pep|seq)$/;
    my $stem = $1;
    open(JOB, ">$jobdir/$stem.pfamscan.sh") || die $!;
    print JOB "#PBS -l nodes=1:ppn=$cpus\n";
    print JOB "$exe --cpu $cpus --tblout $odir/$stem.tbl --domtblout $odir/$stem.domtbl.tab $db $dir/$file > $odir/$stem.$version.hmmscan\n";
    close(JOB);
}
