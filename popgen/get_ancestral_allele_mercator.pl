#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# expecting input to just be tab or space deilmited and first 2 columns are chrom and position

my $mercator_aln_dir = 'alignments';
my $slice_exe = 'sliceAlignment';
my ($ref_genome,$target_genome,$input);
GetOptions('r|ref:s' => \$ref_genome,   # what is the name in the map file?
	   'i|input:s' => \$input,
	   'target:s'  => \$target_genome,
	   'm|mercator:s' => \$mercator_aln_dir,
	   'e|exe:s'   => \$slice_exe, # in case it isn't in your path and you to specify
	   );

die "must provide an input with -i or --input\n" unless defined $input;
die "must provide a ref genome name for when processing the Mercator files\n" unless defined $ref_genome;
die "must provide a target genome name for when processing the Mercator files\n" unless defined $target_genome;


print join("\t", qw(CHROM POS ALLELE)), "\n";
open(my $fh => $input) || die $!;
while(<$fh>) {
    next if /^\#/;

    # all we are looking for is the chrom and position, should be the first 2 columns
    # if not we need to define some other formats
    my ($chrom,$pos,@rest) = split;
    my $exe = sprintf("%s %s %s %s %d %d +", 
		      $slice_exe, $mercator_aln_dir, $ref_genome, $chrom, $pos,$pos+1);
    my $allele = "-";
    open(my $ifh => "$exe 2>&1 |") || die $!;
    while(<$ifh>) {
	if(/^Reading/) { 
	    next;
	} elsif( /^Error:/ ) {
	    last;
	} elsif( /^>\Q$target_genome\E/) {
	    ($allele) = <$ifh>;
	    chomp($allele);
	    last;
	}	
    }
    close($ifh);
    print join("\t" ,$chrom, $pos, $allele),"\n";
}
