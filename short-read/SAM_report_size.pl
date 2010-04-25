#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use List::Util qw(sum);
use Bio::DB::Sam;

# this expects BAM files as the ARGV
my %expected_bases = map { $_ => 1 } qw(C A G T);

my $dir = 'size_summary';
my $minsize = 18;
my $maxsize = 30;
my $genome;
GetOptions('d|dir:s' => \$dir,
	   'min:i'     => \$minsize,
	   'max:i'     => \$maxsize,
	   'g|genome:s'=> \$genome);

mkdir($dir) unless -d $dir;

for my $file ( @ARGV ) {
    my $base = $file;
    $base =~ s/\.bam$//;
    my $db = Bio::DB::Sam->new(-bam => $file,
			       -fasta=>$genome);
    my %counts;
    for my $id ( $db->seq_ids ) {
	my $segment = $db->segment($id);
	my $iterator     = $segment->features(-iterator=>1);
	while (my $align = $iterator->next_seq) { 
	    my $seq = $align->query->dna;
	    my $length = $align->query->length;
	    next if $length < $minsize || $length > $maxsize;
	    my $five_base = uc substr($seq,0,1); # get 5' base
	    my $three_base = uc substr($seq,-1,1); # get 3' base
	    next unless $expected_bases{$five_base};
	    $counts{$length}->{5}->{$five_base}++;
	    next unless $expected_bases{$three_base};
	    $counts{$length}->{3}->{$three_base}++;
	}
    }
    
    open(my $rpt => ">$dir/$base.5_summary_sizes" ) || die $!;
    open(my $rptpct => ">$dir/$base.5_summary_sizes.percent" ) || die $!;
    my @lengths = sort { $a <=> $b } keys %counts;
    my @bases   = sort keys %expected_bases;
    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
    print $rptpct join("\t", qw(LENGTH), @bases),"\n";
    for my $c ( @lengths ) {
	print $rpt join("\t", $c, map { $counts{$c}->{5}->{$_} } @bases), "\n";
	my $sum = sum ( map { $counts{$c}->{5}->{$_} } @bases );
	print $rptpct join("\t", $c, map { sprintf("%.2f",100*$counts{$c}->{5}->{$_}/$sum) } @bases),
	"\n";
    }
    close($rpt);
    close($rptpct);
    open($rpt => ">$dir/$base.3_summary_sizes" ) || die $!;
    open($rptpct => ">$dir/$base.3_summary_sizes.percent" ) || die $!;
    @lengths = sort { $a <=> $b } keys %counts;
    @bases   = sort keys %expected_bases;
    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
    print $rptpct join("\t", qw(LENGTH), @bases),"\n";

    for my $c ( @lengths ) {
	print $rpt join("\t", $c, map { $counts{$c}->{3}->{$_} } @bases), "\n";
	my $sum = sum ( map { $counts{$c}->{5}->{$_} } @bases );
	print $rptpct join("\t", $c, map { sprintf("%.2f",100*$counts{$c}->{3}->{$_}/$sum) } @bases),
	"\n";
    }

}

