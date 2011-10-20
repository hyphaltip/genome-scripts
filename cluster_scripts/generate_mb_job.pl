#!/usr/bin/perl  -w
use Getopt::Long;
use File::Spec;
use strict;

my $odir = ".";
GetOptions('j|job:s' => \$odir);

for my $file (@ARGV) {
 my (undef,$dir,$sfile) = File::Spec->splitpath(File::Spec->rel2abs($file));
 if( $sfile =~ /(\S+)\.nex$/ ) {
  my $ofile = $1;
  open(my $fh => ">$odir/$ofile.sh") ||die $!; 
  print $fh "#PBS -N $ofile\n ",
	"cd $dir\n", "mb $sfile\n";
 close($fh);
 }
}
