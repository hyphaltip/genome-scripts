#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use File::Spec;
use Getopt::Long;
use List::Util qw(sum);
use Bio::SeqIO;
use constant GAP => '-';


# Algorithm
# -- process each alignment, write out score info (gap or identity #) 
#    in genomic coordinates
# -- can dump as a detailed per-site tab file or a wig file 

my ($source,$type) = qw(2way_MERC_PECAN substitution);
my $code = 'Nt'; # call this SNP in context of Ntetrasperma
# variant is 0-based
my $aln_file_name = 'output.mfa';
my $alnformat     = 'fasta';
my $info = '2way_snp';
my $outfile = 'neurospora_crassa_OR74A_7';
my $alndir = 'alignments/pecan_alignments';
my $debug = 0;
my $ref_genome_index = 0;
my @pair = qw(n_crassa n_tetrasperma_2508);
GetOptions(
	   'f|alnfile:s'      => \$aln_file_name,
	   'af|format:s'      => \$alnformat,
	   'v|verbose!'       => \$debug,
	   'o|out|output:s'   => \$outfile,
	   'r|refidx|ref:i'   => \$ref_genome_index,
	   'd|desc:s'         => \$info,
	   'a|i|input:s'      => \$alndir,
	   );

unless( -d $alndir || -r $alndir ) {
    die("cannot open alignment dir '$alndir'\n");
}
open(my $fh => ">$outfile.$info.gff3") || die $!;
print $fh "##gff-version 3\n";
my $mapfile = File::Spec->catfile($alndir,'map');
open(my $mapfh => $mapfile ) || die "$mapfile: $!";
my %chroms;
while(<$mapfh>) {
    my ($aln_id, @line) = split;
    my ($chrom,$start,$end) = map { $line[$ref_genome_index + $_] } 0..2;     
    next if $chrom eq 'NA';    


	my $fname = File::Spec->catfile($alndir,$aln_id,$aln_file_name);
	my $alnio = Bio::AlignIO->new(-file => $fname,
				      -format => $alnformat);
	if( my $aln = $alnio->next_aln ) {
	    my $col_count    = $aln->length;
	    
	    # simplified data structure, we're going to use substr
	    # for speed I hope
	    my (@seqs,@ids);
	    my $i = 0;
	    for my $seq ( $aln->each_seq ) {
		push @seqs, $seq->seq;
		push @ids, $seq->display_id; # just to double check
		$seq = undef;
		last if $i++ == 1;
	    }
	    undef $aln;
	    # process each column
	    my $genome_index = $start+1;
	    for( my $col_index = 0; $col_index < $col_count; $col_index++ ) {
		my $ref_base = uc(substr($seqs[$ref_genome_index],$col_index,1));
		if( $ref_base ne GAP ) { # we skip all ref-base gap columns
		    # as these can't be written in ref coords
		    # iterate thru each sequence
		    my @alleles = ($ref_base);
		    my $i = 0;		
		    for my $rest ( @seqs ) {
			if( $i != $ref_genome_index ) {	# that is NOT the ref
			    next if $pair[$i] ne $ids[$i];
			    my $base = uc(substr($seqs[$i],$col_index,1));
			    if( $base ne GAP && 
				$base ne $ref_base ) {
				# identical
				push @alleles, $base;
			    }
			    last; # short-circuit for now
			}
			$i++;			
		    }
		    if( @alleles > 1) {
			print $fh join("\t", 
				   $chrom,
				   $source,
				   $type,
				   $genome_index,
				   $genome_index,
				   '.',
				   '.',
				   '.',
				   sprintf("ID=%s.%s_%s;Allele=%s;Allele=%s;ref_allele=%s",
					   $code,$chrom,$genome_index,
					   @alleles, $alleles[0])),
			"\n";
		    }
		    $genome_index++;
		}
	    }
	} else {
	    warn("No alignments in $fname\n");
	}
}

