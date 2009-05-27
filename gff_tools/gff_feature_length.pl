#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $type;
my @names = qw(Name ID Parent Note);
GetOptions('t|type:s' => \$type,
	   );
my $total;
while(<>) {
    next if /^\#/;
    chomp;
    my @line = split(/\t/,$_);
    next if( defined $type && $line[2] ne $type );
    my $name = $line[-1];
    my %lastcol = map { split(/=/,$_) } split(/;/,pop @line);
    for my $n( @names ) {
	if( exists $lastcol{$n} ) {
	    $name = $lastcol{$n};
	    last;
	}
    }
    print join("\t", $name, abs($line[4]-$line[3])),"\n";
    $total += abs($line[4] - $line[3]);
}

print " Total = $total\n";
