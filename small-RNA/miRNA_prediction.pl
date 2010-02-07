#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::Sam;
use List::Util qw(sum);
my $min_coverage = 20;
my ($window_table,$bam,$fa);

GetOptions(
	   'w|t|table:s'   => \$window_table,
	   'b|bam:s'       => \$bam,
	   'fa|g|genome:s' => \$fa,
	   'm|min:i'       => \$min_coverage,
	   );

$bam ||= shift @ARGV;
$fa  ||= shift @ARGV;

if( ! defined $window_table ) {
    warn("cannot proceed without a window_table -w/-t\n");
    exit;
}

my $db = Bio::DB::Sam->new(-bam => $bam,
			   -fa  => $fa);
open(my $fh => $window_table ) || die $!;
my $header = <$fh>;
my @header = split(/\t/,$header);

my $count = 0;
while(<$fh>) {
    my ($seqid,$start,$end,@window_counts) = split;
    next unless ( sum(@window_counts) > $min_coverage );
    $count++;
    
}
print "$count windows\n";
