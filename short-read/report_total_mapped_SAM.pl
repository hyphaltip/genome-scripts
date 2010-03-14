#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
use List::Util qw(sum);

# provide SAM file(s) and print out a table of the number of 
# reads binned by windows across the genome

my $window_size = 1000;
my $debug = 0;
my $features;
my $fa;
my $lookup;
GetOptions(
	   'w|window:i' => \$window_size,
 	   'v|verbose!' => \$debug,
	   'l|lookup:s'    => \$lookup,
	   );
my (@bams) = @ARGV;

my %lookups;
if( $lookup ) {
    open(my $fh => $lookup) || die $!;
    while(<$fh>) {
	my @r = split;
	$lookups{$r[0]} = $r[1];
    }
}
my (@dbs,@titles);

my $sam = $dbs[0];
my %counts;
print join("\t",qw(SAMPLE MAPPED_READS)), "\n";

for my $bam ( @bams ) {
    my $title;
    if( $lookups{$bam} ) {
	$title = $lookups{$bam};
    } elsif( $bam =~ /(\S+)\.bam/) {
	$title= $1;
    } else { 
	$title = $bam;
    }

    my $db =  Bio::DB::Sam->new(-bam  =>$bam,
				-fasta=>$fa);

    my $count = 0;
    my @targets    = $db->seq_ids;
    for my $target ( @targets ) {
	my $segment      = $db->segment(-seq_id  => $target);
	my $iterator     = $segment->features(-iterator=>1);
	while( my $aln = $iterator->next_seq ) {
	    if( $aln->name =~ /^(\d+)\-(\d+)$/) {
		$count += $2;
	    } else {
		$count++;
	    }	    
	}
    }
    print join("\t", $title, $count), "\n";
}
