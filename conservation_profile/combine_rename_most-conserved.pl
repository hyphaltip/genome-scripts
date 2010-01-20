#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;

my $ext = 'phastCons_nc_trained_EM.most-conserved.bed';
my $dir = 'phastcons_output';

GetOptions(
	   'e|ext:s' => \$ext,
	   'd|dir:s' => \$dir,
	   );
    
opendir(DIR, $dir) || die $!;
for my $d ( sort { $a->[1] <=> $b->[1] } 
	       map { [$_,split(/\./,$_,2)] } 
	       grep { /^(\d+)\.\Q$ext\E$/ } readdir(DIR) ) {
    my ($file,$num) = @$d;
    open(my $fh => File::Spec->catfile($dir,$file)) || die $!;
    while(<$fh>) {
	s/(output|mavid)/$num/g;
	print;
    }
}
