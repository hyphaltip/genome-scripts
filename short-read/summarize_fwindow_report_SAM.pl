#!/usr/bin/perl -w
use strict;

use Bio::DB::Sam;
use Getopt::Long;

my $debug = 0;
my $features;
GetOptions(
	   'v|verbose!' => \$debug,
	   'f|gff|feature:s' => \$features,
	   );
my ($bam,$fa) = @ARGV;
# high level API
my $sam = Bio::DB::Sam->new(-bam  =>$bam,
			    -fasta=>$fa,
			    );
my @bases = qw(C A G T);


open(my $fh => $features ) || die "$features: $!";

print join("\t", qw(FEATURE TYPE CHROM START END STRAND FWD_STRAND REV_STRAND)), "\n";

while(<$fh>) {
    next if /^\#/;
    my %counts;
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

    my @alignments   = $sam->get_features_by_location(-seq_id => $seqid,
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
	next if keys %f <= 2;	# drop those AAA or TTT runs
	my $skip_this = 0;
	
	my $comp_strand = ($strand == $qstrand) ? 'same' : 'opp';
	
	if( &overlaps($qstart,$qend,$start,$end) ) {
	    my $base = substr($dna,0,1);
	    if( $aln->name =~ /^(\d+)\-(\d+)$/) {
		$counts{$comp_strand} += $2;
	    } else {
		$counts{$comp_strand}++;
	    }
	    last if $debug;
	}	
	last if $debug;
    }
    print join("\t", $name, $type, $seqid, $start,$end, $strand, $counts{'same'} || 0, $counts{'opp'} || 0), "\n";
}


sub overlaps {
    my ($f1_s,$f1_e,
	$f2_s,$f2_e) = @_;
    return not ( $f1_s > $f2_e ||
		 $f1_e < $f2_s);
}
