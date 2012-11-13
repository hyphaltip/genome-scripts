#!/usr/bin/perl -w
use strict;

use Getopt::Long;

my $change_col = 'log2(fold_change)';
my $samp_col_1 = 'sample_1';
my $samp_col_2 = 'sample_2';
my $test_col   = 'significant';
my $in;
GetOptions(
    'c|column:s' => \$change_col,
    's1:s'       => \$samp_col_1,
    's2:s'       => \$samp_col_2,
    'i|in:s'     => \$in,
    'test:s'     => \$test_col,
    );

$in ||= shift || die "need an in file: $!";

open(my $infh => $in) || die $!;
my $header = <$infh>;
chomp($header);
my $i = 0;
my %head = map { $_ => $i++ } split(/\t/,$header);
for my $c ( $change_col, $test_col, $samp_col_1, $samp_col_2) {
    if( ! exists $head{$c} ) {
	die "cannot find expected column $c in header\n";
    }
}
my %data;
my $change_col_n = $head{$change_col};
my $test_col_n   = $head{$test_col};
while(<$infh>) {
    chomp;
    my @r = split(/\t/,$_);
    next unless $r[ $test_col_n ] eq 'yes';
    my ($s1,$s2) = map { $r[ $head{$_} ] } ($samp_col_1,$samp_col_2);
    my $direction = ( $r[ $change_col_n ] < 0 ) ? 'DOWN' : 'UP';
    push @{$data{$direction}->{$s1.'-'.$s2}}, [@r];
}

# sort greatest to least
# use the column name to determine which column to sort by
while ( my ($dir,$d) = each %data ) {
    while( my ($pair,$rows) = each %$d ) {
	open( my $fh => ">$in.$pair\_$dir.dat") || die $!;
	print $fh $header, "\n";
	for my $row ( sort { $b->[ $change_col_n ] <=> $a->[ $change_col_n ] }
		    @$rows ) {
	    print $fh join("\t", @$row), "\n";
	}
    }
}

