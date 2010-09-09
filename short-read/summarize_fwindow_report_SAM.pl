#!/usr/bin/perl -w
use strict;

use Bio::DB::Sam;
use Getopt::Long;

use constant MILLION => 1_000_000;

my $debug = 0;
my $output = "all.mRNA_windows.tab";
my ($fa,$features);
GetOptions(
	   'v|verbose!'       => \$debug,
	   'f|gff|feature:s'  => \$features,
	   'fa:s'             => \$fa,
	   'o|output:s'       => \$output,
	   ); 
my (@bam) = @ARGV;


open(my $allfh => ">$output") || die $!;


# high level API
my %bamfiles;
my @order_bam;
my %perbamfh;
for my $bam ( @bam ) {
    my ($bamname) = $bam;
    $bamname =~ s/\.bam$//;
    push @order_bam, $bamname;
    $bamfiles{$bamname}->{'sam'} = Bio::DB::Sam->new(-bam    => $bam,
						     -fasta  => $fa,
						     );
    
    open(my $perfh => ">$bamname.mRNA_windows.tab") || die $!;
    $perbamfh{$bamname} = $perfh;
    print $perfh join("\t", qw(FEATURE TYPE CHROM START END STRAND FWD REV)), "\n";

    # summary stats for the library for normalization purposes
    open(STATS, "samtools flagstat $bam |") || die $!;
    while(<STATS>) {
	if( /^(\d+)\s+in total/ ) {
	    $bamfiles{$bamname}->{total_reads} = $1;
	} elsif( /^(\d+)\s+mapped/) {
	    $bamfiles{$bamname}->{total_mapped} = $1;
	}	
    }
    printf $allfh "#Library: %s Total: %d Mapped: %d\n", $bamname, $bamfiles{$bamname}->{total_reads}, $bamfiles{$bamname}->{total_mapped};
    $bamfiles{$bamname}->{total_mapped} /= MILLION;
}

my @bases = qw(C A G T);


open(my $fh => $features ) || die "$features: $!";

print $allfh join("\t", qw(FEATURE TYPE CHROM START END STRAND), 
	   map { ($_. '_SAME_',
		  $_. '_OPP_',
		  $_. '_TOTAL_NORM') } @order_bam), "\n";

my $size = 0;

while(<$fh>) {
    next if /^\#/;
    my ($seqid,$src,$type,$start,$end,undef,$strand,undef,$group) = split;
    $strand = -1 if $strand eq '-';
    $strand = 1  if $strand eq '+';

    my %group = map { split(/=/,$_) } split(/;/,$group);
    my $name;
    for my $n ( qw(Id ID Name Parent) ) {
	if( exists $group{$n} ) {
	    $name = $group{$n};
	    last;
	}
    }
    unless ( defined $name ) {
	die("need a parseable group from this line $_\n");
    }
    my @counts;
    for my $bamfile ( @order_bam ) {
	my %count_lib;
	my @alignments   = $bamfiles{$bamfile}->{'sam'}->get_features_by_location(-seq_id => $seqid,
										  -start  => $start,
										  -end    => $end);
	
	for my $aln ( @alignments ) {
	    my $qstart  = $aln->start;
	    my $qend    = $aln->end;
	    my $qstrand = $aln->strand;

	    my $len = $aln->length;
	    warn join("\t", $seqid, $aln->name, $len, $strand),"\n" if $debug;
	    my $dna = $aln->query->dna;
	    my %f;
	    for ( split('',$dna) ) {
		$f{$_}++;
	    }
	    next if keys %f <= 2; # drop those AAA or TTT runs
	    my $skip_this = 0;

	    my $comp_strand = ($strand == $qstrand) ? 'same' : 'opp';

	    if( &overlaps($qstart,$qend,$start,$end) ) {
		my $base = substr($dna,0,1);
		if( $aln->name =~ /^(\d+)\-(\d+)$/) {
		    $count_lib{$comp_strand} += $2;
		} else {
		    $count_lib{$comp_strand}++;
		}
		$count_lib{'total'}++;
		last if $debug;
	    }	
	    last if $debug;	   
	}
	# these are the accumulated counts of number of reads from Same and Opposite strands followed by the number
	# of reads total normalized by the number mapped to the genome
	push @counts, $count_lib{'same'} || 0, $count_lib{'opp'} || 0, sprintf("%.1f",($count_lib{'total'} || 0 ) / 
									       $bamfiles{$bamfile}->{total_mapped});
	my $ofh = $perbamfh{$bamfile};
	print $ofh ("\t", $name, $type, $seqid, $start,$end, $strand, 
				   $count_lib{'same'} || 0, 
				   $count_lib{'opp'} || 0,
				   ), "\n";
    }
    print $allfh join("\t", $name, $type, $seqid, $start,$end, $strand, @counts), "\n";
}


sub overlaps {
    my ($f1_s,$f1_e,
	$f2_s,$f2_e) = @_;
    return not ( $f1_s > $f2_e ||
		 $f1_e < $f2_s);
}
