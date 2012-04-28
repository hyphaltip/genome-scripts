#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
use File::Spec;
use Bio::Perl qw(revcom);
use List::Util qw(sum);

# this expects SOAP output
my %expected_bases = map { $_ => 1 } qw(C A G T);
my %compress = ('bz2' => 'bzcat',
		'gz'  => 'zcat',
		''    => 'cat',
		'Z'   => 'zcat');

my $fa;
my $dir = 'SAM_size_and_partition';
GetOptions('d|dir:s' => \$dir,
	   'fa|g|genome:s' => \$fa,);

my %lookups;
mkdir($dir) unless -d $dir;
for my $bam ( @ARGV ) {
    my $base;
    if( $lookups{$bam} ) {
	$base = $lookups{$bam};
    } elsif( $bam =~ /(\S+)\.bam/) {
	$base = $1;
    } else { 
	$base = $bam;
    }

    my $db =  Bio::DB::Sam->new(-bam  =>$bam,
				-fasta=>$fa,
				);
    my (%ofh,%counts);
    my @targets = $db->seq_ids;
    for my $seq ( map { $db->segment(-seq_id => $_) } @targets ) {
	print $seq->seq_id, " ", $seq->length, " ", $seq->start, "..",$seq->end,"\n";
        for my $read ( $seq->features ) {
	    my $length = $read->query->length;
	    my $ofh = $ofh{$length};
	    unless( defined $ofh ) {
		open($ofh => ">$dir/$base.$length.dat") || die "$dir/$base.$length.dat: $!";
		$ofh{$length} = $ofh;
	    }
	    print $ofh join("\t", $read->seq_id,$read->query->seq_id,
			    $read->start, $read->end, $read->strand,
			    $read->query->start, $read->query->end,
			    $read->query->dna,$read->cigar_str), "\n";
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

