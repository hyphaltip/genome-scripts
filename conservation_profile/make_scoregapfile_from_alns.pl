#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use File::Spec;
use Getopt::Long;
use List::Util qw(sum);
use Bio::SeqIO;
use constant GAP => '-';

my @distance;
$distance[0] = 0; # n.crassa -> n.crassa
$distance[1] = 1; # n.crassa -> n.tetrasperma (=2/3)
$distance[2] = 2; # n.crassa -> n.discreta    (=1/3)
my $Factor = sum(@distance);
for my $x ( @distance ) {
    $x = 1 - $x/$Factor;
}


# Algorithm
# -- process each alignment, write out score info (gap or identity #) 
#    in genomic coordinates
# -- can dump as a detailed per-site tab file or a wig file 

# mercator is 0-based
my $chrom_file    = '/data/genomes/dna/neurospora_crassa_OR74A_7.fasta';
my $aln_file_name = 'output.mfa';
my $alnformat     = 'fasta';
my $info = '3way';
my $outfile = 'wig/neurospora_crassa_OR74A_7';
my $alndir = 'alignments/pecan_alignments';
my $debug = 0;
my $tab = 0; # write out report as TAB or WIG
my $ref_genome_index = 0;

GetOptions(
	   'g|c|genome|chrom:s'=> \$chrom_file,
	   'f|alnfile:s'      => \$aln_file_name,
	   'af|format:s'      => \$alnformat,
	   'v|verbose!'       => \$debug,
	   't|tab'            => \$tab, # TAB or WIG flag
	   'o|out|output:s'   => \$outfile,
	   'r|refidx|ref:i'   => \$ref_genome_index,
	   'd|desc:s'         => \$info,
	   'a|i|input:s'      => \$alndir,
	   );

my $gapname = "$info\_gaps";
my $gapdesc = "$info gapped alignment calc";

my $idname = "$info\_identical";
my $iddesc = "$info identical alignment calc";

unless( -d $alndir || -r $alndir ) {
    die("cannot open alignment dir '$alndir'\n");
}
my (%chrom_lengths,@chrom_order);
{
    my $inseq = Bio::SeqIO->new(-format => 'fasta',
				-file   => $chrom_file);
    while( my $seq = $inseq->next_seq ) {
	$chrom_lengths{$seq->display_id} = $seq->length;
	push @chrom_order, $seq->display_id;
    }
}
my $mapfile = File::Spec->catfile($alndir,'map');
open(my $mapfh => $mapfile ) || die "$mapfile: $!";
my %chroms;
while(<$mapfh>) {
    my ($aln_id, @line) = split;
    my ($chrom,$start,$end) = map { $line[$ref_genome_index + $_] } 0..2;     
    next if $chrom eq 'NA';    
    push @{$chroms{$chrom}}, [$aln_id, $start,$end]
} 

my ($gapfh,$idfh);
if( $tab ) {
    open($gapfh => ">$outfile.$info.tab") || die $!;
    print $gapfh join("\t", qw(CHROM POSITION 
			       ALN_ID ALN_COL 
			       ALLELES IDENTICAL GAPPED)),"\n";
} else {
    open($gapfh => ">$outfile.$info\_gap.wig") || die $!;
    open($idfh => ">$outfile.$info\_identity.wig") || die $!;
    printf $gapfh "track type=wiggle_0 name=\"%s\" description=\"%s\"\n",
    $gapname,$gapdesc;
    printf $idfh "track type=wiggle_0 name=\"%s\" description=\"%s\"\n",
    $idname,$iddesc;
}

for my $chrom ( @chrom_order ) { 
    my (%scores,%gaps,%info);     
    unless( $tab ) {
	printf $gapfh "fixedStep chrom=%s start=%s step=1\n",$chrom,1;
	printf $idfh  "fixedStep chrom=%s start=%s step=1\n",$chrom,1;
    }
    for my $pos ( @{$chroms{$chrom}} ) {
	my ($aln_id, $start,$end) = @$pos;

	my $fname = File::Spec->catfile($alndir,$aln_id,$aln_file_name);
	my $alnio = Bio::AlignIO->new(-file => $fname,
				      -format => $alnformat);
	if( my $aln = $alnio->next_aln ) {
	    my $col_count    = $aln->length;

	    # simplified data structure, we're going to use substr
	    # for speed I hope
	    my (@seqs,@ids);
	    for my $seq ( $aln->each_seq ) {
		push @seqs, $seq->seq;
		push @ids, $seq->display_id; # just to double check
		$seq = undef;
	    }
	    undef $aln;
	    # process each column
	    my $genome_index = $start+1;
	    for( my $col_index = 0; $col_index < $col_count; $col_index++ ) {
		my $ref_base = uc(substr($seqs[$ref_genome_index],$col_index,1));
		if( $ref_base ne GAP ) { # we skip all ref-base gap columns
		    # as these can't be written in ref coords
		    # iterate thru each sequence
		    my @debug = ($ref_base);
		    my $i = 0;		
		    my ($score,$gap) = (0,0);
		    for my $rest ( @seqs ) {
			if( $i != $ref_genome_index ) {	# that is NOT the ref
			    my $base = uc(substr($seqs[$i],$col_index,1));
			    if( $base eq $ref_base ) {
				# identical
				$score += $distance[$i];
			    } elsif( $base eq GAP ) {
				$gap   += $distance[$i];
			    }
			    push @debug, $base;
			}
			$i++;			
		    }
		    $scores{$genome_index} = $score;
		    $gaps{$genome_index}   = $gap;
		    if( $tab ) {
			$info{$genome_index} = [ $aln_id,
						 $col_index,
						 join(",",@debug),
						 ];
		    }
		    $genome_index++;
		}
	    }
	} else {
	    warn("No alignments in $fname\n");
	}
	warn("  done with $chrom ($aln_id)\n");
    }
    warn("done with $chrom\n"); 
    for( my $i = 1; $i <= $chrom_lengths{$chrom}; $i++ ) {
	if( $tab ) {	    
	    if( exists $info{$i} ) {
		print $gapfh join("\t", $chrom, $i, @{$info{$i}},
				  $scores{$i},$gaps{$i}),"\n";
	    } 
	} else {
	    if( exists $gaps{$i} ) {
		printf $gapfh "%.2f\n", 1 - $gaps{$i};
		printf $idfh "%.2f\n", $scores{$i};
	    } else {
		print $gapfh "1.00\n";
		print $idfh "0.00\n";
	    }
	}
    }
}

