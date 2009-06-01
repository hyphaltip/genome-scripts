#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max);
my $only_type;

GetOptions('t|type:s' => \$only_type);
my ($category, $data) = @ARGV;

open(my $fh => $category) || die $!;
my %cat_data;
while(<$fh>) {
    next if /^\#/;
    chomp;
    my ($seqid,$src,$type,$start,$end,$score,$strand,$frame,
	$group) = split(/\t/,$_);
    next if( defined $only_type && $type ne $only_type);

    my %grp = map { split(/=/, $_) } split(/;/,$group);
    my $id = $grp{'ID'};
    if( $grp{'Name'} ) {
	$id .= "_".$grp{'Name'};
    }
    push @{$cat_data{$seqid}}, [(sort { $a <=> $b } $start,$end),$id];
}

my %counts;
open($fh => $data) ||die $!;
my $i;
while(<$fh>) {
    next if /^\#/;
    chomp;
    my ($seqid,$src,$type,$start,$end,$score,$strand,$frame,
	$group) = split(/\t/,$_);
    if( exists $cat_data{$seqid} ) {
	for my $c ( @{$cat_data{$seqid}} ) {
	    if( overlaps([sort { $a <=> $b } $start,$end], 
			 [$c->[0],$c->[1]]) ) {
		$counts{$seqid}->{$c->[2]}++;
	    }
	}
    }   
}
for my $s ( sort keys %counts ) {
    print "$s\n";
    for my $r ( sort keys %{$counts{$s}} ) {
	print "\t$r\t",$counts{$s}->{$r},"\n";
    }
}
sub overlaps {
    my ( $left, $right ) = @_;
    return ! ( $left->[0] > $right->[1] ||
	       $left->[1] < $right->[0]);
}

