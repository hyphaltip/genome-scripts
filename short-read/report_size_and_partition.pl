#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
my %expected_bases = map { $_ => 1 } qw(C A G T);
my %compress = ('bz2' => 'bzcat',
		'gz'  => 'zcat',
		''    => 'cat',
		'Z'   => 'zcat');

my $dir = 'size_summary';
GetOptions('d|dir:s' => \$dir);

mkdir($dir) unless -d $dir;

for my $file ( @ARGV ) {
    my (%ofh,%counts);
    my (undef,undef,$fname) = File::Spec->splitpath($file);
    my ($base,$ext) = split(/\.([^\.]+)$/,$fname,2);
    my $fh;
    if( $ext ) {
	open($fh => "$compress{$ext} $file |") || die $!;
    } else {
	open($fh => "<$file") || die $!;
    }
    while(<$fh>) {
	chomp;
	my ($name,$seq,$qual,$hitcount,
	    $flag, $length, $direction,
	    $chrom, $location_start, # 1 based
	    $types) = split(/\t/,$_);

	my $ofh = $ofh{$length};
	unless( defined $ofh ) {
	    open($ofh => ">$dir/$base.$length.dat") || die $!;
	    $ofh{$length} = $ofh;
	}
	print $ofh $_,"\n";
	my $five_base = uc substr($seq,0,1); # get 5' base
	next unless $expected_bases{$five_base};
	$counts{$length}->{$five_base}++;
    }

    open(my $rpt => ">$dir/$base.summary_sizes" ) || die $!;
    my @lengths = sort { $a <=> $b } keys %counts;
    my @bases   = sort keys %expected_bases;
    print $rpt join("\t", qw(LENGTH), @bases),"\n";
    for my $c ( @lengths ) {
	print $rpt join("\t", $c, map { $counts{$c}->{$_} } @bases), "\n";
    }
}
