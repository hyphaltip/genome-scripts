#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Getopt::Long;
use Bio::AlignIO;

my $slicer = 'sliceAlignment';
my $target = 'Ccin';
my $window = 120;
my $overlap = 40;
my $idir = 'alignments';
my $odir = 'output_slices';
my $jobdir = 'jobs';

GetOptions(
	   'e|exe:s'        => \$slicer,
	   'i|in|input:s'   => \$idir,
	   'o|out|output:s' => \$odir,
	   'j|job|jobdir:s' => \$jobdir,
	   );


for my $d ( $jobdir,$odir) {
    mkdir($d);
}
open(my $genomes => "$idir/genomes") || die $!;
my %genome2col;
while(<$genomes>) {
    my $i =0;
    %genome2col = map { $_ => $i++ } split;
    last;
}
my $target_col;
if( exists $genome2col{$target} ) {
    $target_col = $genome2col{$target};
} else {
    die("Target $target is unknown genome, examine $idir/genomes for the expected names\n");
}
open(my $fh => "$idir/map") || die $!;
my @intervals;
while(<$fh>) {
    my ($int,@data) = split;
    push @intervals, [ splice(@data,$target_col*4,4)];
}
for my $interval ( @intervals ) {
    my ($chrom,$begin,$end,$strand) = @$interval;
    open(my $jobfh => ">$jobdir/$chrom.$begin-$end.sh") || die $!;
    for( my $start = $begin; $start < $end; $start+= $overlap) {
	my $stop = $start + $window;
	my $interval_name = sprintf("%s_%d-%d",$chrom,$start,$stop);
	open(my $ifh => "sliceAlignment $idir $target $chrom $start $stop + |") || die;
	my $alnIO = Bio::AlignIO->new(-fh => $ifh,
				      -alphabet=>'dna',
				      -format => 'fasta');
	my $aln = $alnIO->next_aln;
	if( $aln && $aln->length > 0 ) {
	    my $min_len = 70;
	    my $keep = 1;
	    my $aln_len = $aln->length;
	    for my $seq ( $aln->each_seq ) {		
		$keep=0,last if( ($aln_len - ($seq->seq=~ tr/-/-/)) 
				 < $min_len);
	    }
	    if( $keep ) { 
		$aln->set_displayname_flat(1);
		my $alnout = Bio::AlignIO->new(-format => 'clustalw',
					       -alphabet=>'dna',
					       -file   => ">$odir/$interval_name.aln");
		$alnout->write_aln($aln);
		print $jobfh "RNAz $odir/$interval_name.aln > $odir/$interval_name.RNAz\n";
	    }
	}
    }
    close($jobfh);
}
