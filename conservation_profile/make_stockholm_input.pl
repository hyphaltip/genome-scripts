#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Bio::AlignIO;
use Bio::TreeIO;

my $debug   = 0;
my $format  = 'fasta';
my $aln_name= 'mavid.mfa'; # hardcoded to mavid.mfa for now
my $dir     = 'alignments/pecan_alignments';
my $odir    = 'input';
my $ext     = 'pecan';
GetOptions(
           'v|verbose|debug!' => \$debug,
	   'f|format:s' => \$format,
	   'n|name:s'   => \$aln_name,
	   'd|dir:s'    => \$dir,
	   'o|odir:s'   => \$odir,
	   'e|ext:s'    => \$ext,
	   );

mkdir($odir) unless -d $odir;

opendir(my $dirfh => $dir) || die "$dir: $!";
for my $d ( readdir($dirfh) ) {
    next if $d =~ /^\./;
    next unless -d "$dir/$d";
    my $treefh;
    unless( open($treefh => "$dir/$d/treefile") ) {
	warn "$dir/$d/treefile: $!";
	next;
    }
    my $treestring = <$treefh>;
	# fix trailing label on root node that is somehow broken
   $treestring =~ s/:0;$/;/;
    my $aln = Bio::AlignIO->new(-format => $format,
				-file   => "$dir/$d/$aln_name")->next_aln;
    
    open(my $fh => ">$odir/$d.$ext.stock") || die $!;
    print  $fh "# STOCKHOLM 1.0\n";
    printf $fh "%-10s %s\n",'#=GF NH',$treestring;
    for my $seq ( $aln->each_seq ) {
	printf $fh "%-10s %s\n", $seq->display_name, uc($seq->seq);
    }
    print $fh "//\n";
    last if $debug;
}
