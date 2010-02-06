#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;
use Env qw(USER HOME);
use Getopt::Long;
use File::Spec;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use Data::Dumper;

my $debug;
my $dir = 'promotor_alignments';
my @sp  = qw(Ncra Ntet Ndis);
my $genome = 'Smac';
GetOptions(
	   'i|in|dir:s' => \$dir,
	   'g|genome:s' => \$genome,
	   'v|verbose!'=> \$debug,
	   );

my $stats = Bio::Align::DNAStatistics->new;

opendir(DIR, $dir) || die $!;

for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.fasta$/;
    my $stem = $1;
    my $alndir = File::Spec->catfile($dir,"$stem.d");
    opendir(ALNDIR, $alndir) || die $!;
    open(my $fh => ">$dir/$stem.stats.dat") || die $!;
    print $fh join("\t", qw(GENE ALL_ID NC_ID), map { uc($_ ."_D") } @sp),"\n";
    for my $fd ( readdir(ALNDIR) ) {
	next unless ( $fd =~ /(\S+)\.fas$/);
	my $fstem = $1;
	my $in = File::Spec->catfile($alndir,$fd);
	my $input = Bio::AlignIO->new(-format => 'fasta',
				      -file   => $in);
	my $aln = $input->next_aln;
	next unless $aln;
	my $dist = $stats->distance(-align => $aln,
				    -method => 'k80');
	next unless( defined $dist );

	my $uncorr_dist = $stats->distance(-align => $aln,
					   -method => 'uncorrected');
	print $fh join("\t", $fstem, $aln->average_percentage_identity,
		       $dist->get_entry($genome,'Ncra'),
		       map { $dist->get_entry($genome,$_) } @sp),"\n";
		   
	last if $debug;
    }
}
