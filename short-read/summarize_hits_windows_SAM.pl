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
my $lookup; # a summary table of short names to long names of strain files
GetOptions( 
	   'w|window:i' => \$window_size,
	   'v|verbose!' => \$debug,
	   'f|gff|feature:s' => \$features,
	   'fa|g|genome:s' => \$fa,
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
for my $bam ( @bams ) {
    if( $lookups{$bam} ) {
	push @titles, $lookups{$bam};
    } elsif( $bam =~ /(\S+)\.bam/) {
	push @titles, $1;
    } else { 
	push @titles, $bam;
    }

    push @dbs,  Bio::DB::Sam->new(-bam  =>$bam,
				  -fasta=>$fa,
				  );
}
my %skip;
if( $features ) {
    open(my $fh => $features ) || die $!;
    while(<$fh>) {
	next if /^\#/;
	my ($seqid,$src,$type,$start,$end) = split;
	push @{$skip{$seqid}}, [ $start, $end];
    }
}
my $sam = $dbs[0];
my %counts;
my @targets    = $sam->seq_ids;
print join("\t", qw(SEQID START END), @titles), "\n";

for my $target ( @targets ) {
    my $len = $sam->length($target);
    for( my $window = 1; $window < $len; $window += $window_size ) {
	my $end = $window+$window_size-1;
	if( $end > $len ) {
	    $end = $len;
	}
	my @vals;
	for my $db ( @dbs ) {
	    my $count = 0;
	    for my $aln ($db->get_features_by_location
		       (-seq_id => $target,
			-start  => $window,
			-end    => $end) ) {
		if( $aln->name =~ /^(\d+)\-(\d+)$/) {
		    $count += $2;
		} else {
		    $count++;
		}	    
	    }
	    push @vals, $count;
	}
	next if sum( @vals) == 0;
	print join("\t", $target, $window,$end,@vals), "\n";
    }
}
