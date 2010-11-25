#!/usr/bin/perl -w
use strict;

use Bio::DB::Sam;
use Getopt::Long;

my %expected_bases = map { $_ => 1 } qw(C A G T);
my @bases   = sort keys %expected_bases;
my @strands = qw(same opp);

my ($minsize,$maxsize) = (17,30);
my $debug = 0;
my $features;
my $genome;
my $odir = 'feature_read_count';
GetOptions(
	   'min:i'           => \$minsize,
	   'max:i'           => \$maxsize,
	   'd|dir:s'         => \$odir,
	   'v|verbose!'      => \$debug,
	   'f|gff|feature:s' => \$features,
	   'g|fa|genome:s'   => \$genome,
	   );
die("must have feature file") unless defined $features && -f $features;
die("must have genome fasta") unless defined $genome && -f $genome;

mkdir($odir) unless -d $odir;
open(my $fh => $features ) || die "$features: $!";

my (%collected,%totals);
my @bams = map { [$_ =~ /^([^\.]+)\./, Bio::DB::Sam->new(-bam => $_, -fasta => $genome)] } @ARGV;

while(<$fh>) {
    next if /^\#/;
    my ($seqid,$src,$type,$start,$end,undef,$strand,undef,$group) = split;
    $strand = -1 if $strand eq '-';
    $strand = 1  if $strand eq '+';

    my %group = map { split(/=/,$_) } split(/;/,$group);
    my $fname;
    for my $n ( qw(Id ID Name Parent) ) {
	if( exists $group{$n} ) {
	    $fname = $group{$n};
	    last;
	}
    }
    unless ( defined $fname ) {
	die("need a parseable group from this line $_\n");
    }

    for my $bamd ( @bams ) {
	my ($dbname,$db) = @$bamd;
	my $segment = $db->segment($seqid, $start => $end);

	my $iterator = $segment->features(-iterator => 1);
    	while (my $aln = $iterator->next_seq) {				    
	    my $dna = $aln->query->dna;
	    my %f;
	    for ( split('',$dna) ) { $f{$_}++ }
	    next if keys %f <= 2; # drop those AAA or TTT runs
	    my $n = 1;
	    if( $aln->name =~ /^(\d+)\-(\d+)$/) {
		$n = $2;
	    }
#	    warn("seq is ",$aln->query->seq_id," ", $dna,"\n");

	    my $qstrand = $aln->query->strand;
	    my $comp_strand = ($strand == $qstrand) ? 'same' : 'opp';
	    my $five_base = uc substr($dna,0,1);
	    next unless $expected_bases{$five_base};
	    my $len = $aln->length;
	    next if $len < $minsize || $len > $maxsize;
	    my $three_base = uc substr($dna,-1,1);
	    next unless $expected_bases{$three_base};
#	    warn("  n is $n\n");
	    $collected{$dbname}->{$fname}->{$len}->{'all'}->{$comp_strand} += $n;
	    $collected{$dbname}->{$fname}->{$len}->{'5'}->{$five_base} += $n;   
	    $collected{$dbname}->{$fname}->{$len}->{'3'}->{$three_base} += $n;    	    
	    $totals{$dbname}->{$fname} += $n;
	}		
	last if $debug;
    }
}
close($fh);
for my $bamdb ( keys %totals ) {
    open(my $ofh => ">$odir/$bamdb.counts.dat") || die $!;	
    print $ofh join("\t",qw(FEATURE LENGTH), (map { uc $_ }  @strands), 
		    (map { $_."_5" } @bases),
		    (map { $_ . "_3" } @bases)), "\n"; 	
    for my $feature ( sort { $totals{$b} <=> $totals{$a} } keys %{$totals{$bamdb}}){
	for my $len ( sort { $a <=> $b } keys %{$collected{$bamdb}->{$feature}} ) {
	    print $ofh join("\t", $feature, $len,
			    (map { $collected{$bamdb}->{$feature}->{$len}->{all}->{$_} || 0} @strands),
			    (map { $collected{$bamdb}->{$feature}->{$len}->{5}->{$_} || 0} @bases ),
			    (map { $collected{$bamdb}->{$feature}->{$len}->{3}->{$_} || 0 } @bases )), "\n";
	}
    }
}

