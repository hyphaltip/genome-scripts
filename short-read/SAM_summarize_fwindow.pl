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

my %targets;
open(my $fh => $features ) || die "$features: $!";

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
	warn join("\t", $seqid, $aln->name, $len, $strand),"\n";#if $debug;
	my $dna = $aln->query->dna ;
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
		$counts{$comp_strand}->{$len}->{$base} += $2;
		$total+=$2;
	    } else {
		$counts{$comp_strand}->{$len}->{$strand}->{$base}++;
		$total++;
	    }
	    last if $debug;
	}	
	last if $debug;
    }
    warn("total is $total\n");
    for my $dirw ( keys %counts ) {
	open(my $ofh => ">$name.$dirw.tab") || die $!;	
	print $ofh join("\t",qw(LENGTH), @bases), "\n"; 
	for my $size ( sort { $a <=> $b } keys %{$counts{$dirw}} ) {
	    print $ofh join("\t", $size, map { $counts{$dirw}->{$size}->{$_} || 0 } @bases),"\n";
	}
    }
}


sub overlaps {
    my ($f1_s,$f1_e,
	$f2_s,$f2_e) = @_;
    return not ( $f1_s > $f2_e ||
		 $f1_e < $f2_s);
}
