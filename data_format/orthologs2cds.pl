#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Spec;
use Bio::DB::Fasta;

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat');

my $seqdir = 'cds';
my $type   = 'cds';
my $outdir = 'orthologs';
my $prefix = 'orthmcl_';
my $in;
my $single_copy = 1;
GetOptions(
	   'singlecopy!'=> \$single_copy,
	   's|seqdir:s' => \$seqdir,
	   'o|outdir:s' => \$outdir,
	   'p|prefix:s' => \$prefix,
	   't|type:s'   => \$type,
	   'i|input:s'  => \$in,
	   );

my $infh;
if( $in ) {
    if( $in =~ /(?:\.(gz|bz2))$/ ) {
	open($infh, "$uncompress{$1} $in |") || die "cannot open $uncompress{$1} $in:$!";
    } else {
	open($infh, "< $in") || die "cannot open $in:$!";
    }
} else {
    $infh = \*ARGV;
}
mkdir($outdir) unless -d $outdir;
my $db = Bio::DB::Fasta->new($seqdir);
while(<$infh>) {
    my ($orthid, @genes);
    if( s/^ORTHOMCL(\d+)\((\d+) genes,(\d+) taxa\):\s+//) { 

	$orthid = $1;
	my ($genecount,$taxacount) = ($2,$3);
	next if( $single_copy && $genecount != $taxacount); # insure single-copy orthologs
	s/\(\S+\)//g;
	@genes = split(/\s+/,$_);
    } else {
	warn;
	($orthid,@genes) = split;
    }

    my $filename = sprintf("%s%d.%s.fa",$prefix,$orthid,$type);
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
				-file   => '>'.
				File::Spec->catfile($outdir,$filename));
    for my $g ( @genes ) {
	my ($id,$desc) = split(/\s+/,$db->header($g),2);
	    
	my $seq = Bio::PrimarySeq->new(-seq  => $db->seq($g),
				       -id   => $g,
				       -desc => $desc);
	$seqio->write_seq($seq);
    }
}



