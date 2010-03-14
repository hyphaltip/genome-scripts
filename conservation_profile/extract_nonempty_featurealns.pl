#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use File::Spec;
use File::Copy qw(move);
use Getopt::Long;

# reject alignments that are missing any species
my $dry_run = 0;
my $dir;
my $skip_dir = 'skip';
my $ext = 'fas';
my $format = 'fasta';
my $char_ratio_cutoff = '0.75'; # every seq must participate 75%
GetOptions(
	   'f|format:s'=> \$format,
	   'd|dir=s' => \$dir,
	   'e|ext:s' => \$ext,
	   'dryrun!' => \$dry_run,
	   's|skip:s'=> \$skip_dir,
	   'c|cutoff:f'=> \$char_ratio_cutoff);
die("need a dir with -d\n") unless $dir && -d $dir;
$skip_dir = File::Spec->catdir($dir,$skip_dir);
mkdir($skip_dir) unless -d $skip_dir;
	   
opendir(my $DIR => $dir) || die $!;
for my $file ( readdir($DIR) ) {
    next unless $file =~ /(\S+)\.\Q$ext\E$/;
    my $infile = File::Spec->catfile($dir,$file);
    my $in = Bio::AlignIO->new(-format => $format,
			       -file   => $infile,
			       -alphabet => 'dna');
    my $skip_file = 0;
    if( my $aln = $in->next_aln ) {
	my $aln_len = $aln->length;
	for my $seq ( $aln->each_seq ) {
	    my $str = $seq->seq;
	    my $gap = $str =~ tr/-/-/;
	    my $char_ratio = 1 - ($gap / $aln_len);
	    if( $char_ratio < $char_ratio_cutoff) {
		$skip_file = 1;
		last; # short circuit if any sequence is < expected ratio
	    }

	}
    }
    if( $skip_file ) {
	my $rename = File::Spec->catfile($skip_dir,$file);
	if( $dry_run ) {
	    warn ("mv $infile $rename\n");
	} else {
	    move($infile,$rename);
	}
    }
}
