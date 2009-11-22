#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($src,$type) = qw(syntney Lbic);

GetOptions(
	   's|src:s' => \$src,
	   't|type:s'=> \$type,
	   );

my %seen;
my @lookup = qw(0 I II III IV V VI VII VIII IX X XI XII XIII);
while(<>)  {
    next if /^\#/;
    my @line = split;
    
    if( $line[0] =~ /Chr_(\d+)/) { 
	$line[0] = sprintf("Chr_%s",$lookup[$1]);
    }
    if( $seen{$line[3]}++ ) {
	$line[3] .= "_".$seen{$line[3]};
    }
    print join("\t", $line[0], $src, $type, $line[1],$line[2], ".","+",".","ID=$line[3]"),"\n";
}
