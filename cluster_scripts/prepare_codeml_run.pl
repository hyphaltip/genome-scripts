#!/usr/bin/perl -w
use strict;
use File::Spec;
my $jobdir = 'jobs';

my $dir = shift || die $!;
$dir = File::Spec->rel2abs($dir);

mkdir($jobdir);
opendir(DIR, $dir) || die $!;
my $CODEML = '/srv/projects/stajichlab/bin/codeml';

my @Site = qw(m0 m1Neutral m2Selection m3Discrtk2 m3Discrtk3 m7 m8 m8a);
my @Branch = qw(modelA modelAnull modelB);
for my $fdir ( readdir(DIR) ) {
    next unless $fdir =~ /(\S+)\.cds$/;
    my $stem = $1;
    open(my $fh =>">$jobdir/$stem.sh");
    print $fh "#!/bin/bash\n";
    for my $s ( @Site ) {
	print $fh "cd $dir/$fdir/SitePAML/$s\n", "$CODEML codeml.ctl >& codeml.out\n";
    }	
for my $s ( @Branch ) {
    print $fh "cd $dir/$fdir/BranchPAML/$s\n", "$CODEML codeml.ctl >& codeml.out\n";
    }
    close($fh);
}
