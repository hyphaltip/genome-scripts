#!/usr/bin/perl -w
# $Id$

=head1 NAME

gff_feature_length - summary of feature lengths in a GFF file

=head1 USAGE

gff_feature_length file.gff3

=head1 DESCRIPTION 

 Will print out per-feature lengths and total length 

CMDLINE options
 -type provides ability to filter only by certain features (i.e. gene, mRNA)

=cut

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
