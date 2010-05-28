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
my $print_name = 0;
my @names = qw(Name Parent ID Note);
GetOptions('t|type:s' => \$type,
	   'g|genome:i'=> \$genome,
	   'p|print|v!'=> \$print_name,
	   );
my $total;
my $count;
my %transcript;
my $stats = Statistics::Descriptive::Full->new;
while(<>) {
    next if /^\#/;
    chomp;
    my @line = split(/\t/,$_);
    my $name = $line[-1];    
    my %lastcol = map { split(/=/,$_) } split(/;/,pop @line);
    
    if( $line[2] eq 'exon' ) {
	if( exists $lastcol{'Parent'} ) {
	    $transcript{$lastcol{'Parent'}} += abs($line[4]-$line[3])+1;
	} else {
	    die("could not determine a group for this exon\n");
	}
    }
    next if( defined $type && $line[2] ne $type );
    for my $n ( @names ) {
	if( exists $lastcol{$n} ) {
	    $name = $lastcol{$n};
	    last;
	}
    }
    my $len = abs($line[4]-$line[3])+1;
    print join("\t", $name, ),"\n" if $print_name;
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
$stats = Statistics::Descriptive::Full->new;
$stats->add_data([values %transcript]);

printf "Transcript data Max = %d Mean = %.1f Median = %1.f Total = %d N = %d\n",
    $stats->max, $stats->mean, $stats->median,$stats->sum,$stats->count;

if( $genome ) {
    printf "  %.2f Mb Total %% %.2f \n",$stats->sum / 1000000, 
    100 * $stats->sum / $genome;
}
