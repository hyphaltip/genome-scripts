#!/usr/bin/perl -w
use strict;

use Getopt::Long;

my $col = 'FPKM';

GetOptions(
    'c|column:s' => \$col,
    );

my $header = <>;
chomp($header);
my $i = 0;
my %head = map { $_ => $i++ } split(/\t/,$header);
my @data;
while(<>) {
    chomp;
    push @data, [split(/\t/,$_)];
}

# sort greatest to least
# use the column name to determine which column to sort by
print $header, "\n";
for my $d ( sort { $b->[ $head{$col} ] <=> $a->[ $head{$col} ] }
	    @data ) {
    print join("\t", @$d), "\n";
}

