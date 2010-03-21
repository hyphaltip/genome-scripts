#!/usr/bin/perl -w
use strict;

use Bio::DB::Fasta;

my $genome =shift;
my $report = shift;


my $dbh = Bio::DB::Fasta->new($genome);
my $genome_len;
my $stream = $dbh->get_PrimarySeq_stream;
while( my $seq = $stream->next_seq ) {
    $genome_len += $seq->length;    
}
my %chroms;
open(my $fh => $report) || die $!;
my $n = 0;
while(<$fh>) {
    next if /^\s+(SW|score)/ || /^\s+$/ ;
    my ($score,$div,$del,$ins,$query,$qstart,$qend,$q_left,$strand,
	$family, $class, $hstart,$hend,$h_left,$id) = split;
    if( ! exists $chroms{$query} ) {
	$chroms{$query}->[$dbh->length($query)] = {};
    }
    for my $i ( $qstart..$qend ) {
	$chroms{$query}->[$i]->{"$family,$class"}++;
    }
    $n++;
}

my %totals_class;
my %totals_family;
while( my ($chrom,$bases) = each %chroms ) {
    for my $b ( @$bases ) {
	next unless defined $b;
	for my $class ( keys %{$b} ) {
	    my ($family,$cl) = split(',',$class);
	    $totals_class{$cl}++;
	    $totals_family{$class}++;
	}
    }
}

print "Families\n";
for my $type ( sort { $totals_family{$b} <=> $totals_family{$a} } 
	       keys %totals_family ) {
    print join("\t", $type, $totals_family{$type}, 
	       sprintf("%.2f",100 * ( $totals_family{$type} / $genome_len))),
    "\n";
}

print "\nClasses\n";
for my $type ( sort { $totals_class{$b} <=> $totals_class{$a} } 
	       keys %totals_class ) {
    print join("\t", $type, $totals_class{$type}, 
	       sprintf("%.2f",100 * ( $totals_class{$type} / $genome_len))),
    "\n";
}
