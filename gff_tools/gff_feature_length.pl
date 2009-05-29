#!/usr/bin/perl -w
# $Id$

=head1 NAME

gff_feature_length - summary of feature lengths in a GFF file

=head1 USAGE

gff_feature_length file.gff3

=head1 DESCRIPTION 

 Will print out per-feature lengths and total length 

CMDLINE options
 -type [type] provides ability to filter only by certain features (i.e. gene, mRNA)
 -genome [size] genome size to give fraction of genome
=cut

use strict;
use Getopt::Long;
use Statistics::Descriptive;

my $type;
my $genome;
my @names = qw(Name ID Parent Note);
GetOptions('t|type:s' => \$type,
	   'g|genome:i'=> \$genome,
	   );
my $total;
my $count;
my $stats = Statistics::Descriptive::Full->new;
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
    my $len = abs($line[4]-$line[3])+1;
    print join("\t", $name, ),"\n";
    $total += $len;
    $stats->add_data($len);
    $count ++;
}

printf " Max = %d Mean = %.1f Median = %1.f Total = %d N = %d\n",
    $stats->max, $stats->mean, $stats->median,$stats->sum,$stats->count;

if( $genome ) {
    printf "  %.2f Mb Total %% %.2f \n",$stats->sum / 1000000, 
    100 * $stats->sum / $genome;
}
