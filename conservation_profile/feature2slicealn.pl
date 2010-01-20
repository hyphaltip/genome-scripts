#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;
use Env qw(USER HOME);
use Getopt::Long;
use File::Spec;
my $alndir = File::Spec->catfile('alignments','pecan_alignments');
my $genome = 'n_crassa';
my $debug;
my $dir = 'feature_alignments';
GetOptions(
	   'i|in|dir:s' => \$dir,
	   'a|align:s'  => \$alndir,
	   'g|genome:s' => \$genome,
	   'v|verbose!'=> \$debug,
	   );

opendir(DIR, $dir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.fasta$/;
    warn("file is $file\n");
    my $stem = $1;
    my $odir = File::Spec->catfile($dir,"$stem.d");
    mkdir($odir);
    my $in = File::Spec->catfile($dir,$file);
    open(my $fh => "grep '^>' $in |") || die $!;
    while(<$fh>) {
	if( /^>(\S+)\s+(\S+)/) {
	    my ($id,$loc) = ($1,$2);
	    my ($chrom,$pos) = split(/:/,$loc);
	    my $strand = '+';
	    if( $loc =~ s/complement\(// ) {
		$strand = '-';
		$loc =~ s/\)$//;
	    }
	    my ($start,$end) = ($loc =~ /(\d+)\.\.(\d+)/);
	    $start--; #
	    `sliceAlignment $alndir $genome $chrom $start $end $strand > $odir/$id.fas`;
	}
	last if $debug;
    }
}
