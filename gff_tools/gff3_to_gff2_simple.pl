#!/usr/bin/perl -w
use strict;

my $group_name = 'GenePrediction';
my %names;

while( <> ) {
    next if /^\#/;
    chomp;
    my @row = split(/\t/,$_);
    my %rest = map { split(/=/,$_) } split(/;/,pop @row);
    
    if( $row[2] eq 'CDS' || $row[2] eq 'cds' ) {
	if( exists $rest{'Parent'} ) {
	    my $group;
	    if( exists $names{$rest{'Parent'}} ){ 
		$group = $names{$rest{'Parent'}};
	    } else {
		$group = $rest{'Parent'};
	    }
	    print join("\t", @row, 
		       sprintf("%s %s",$group_name,$group)),"\n";	  
	} else {
	    warn("Missing required 'Parent' tag\n");
	}
    } elsif( $row[2] eq 'mRNA' ) {
	my $name;
	for my $key ( qw(Name ID Parent) ) {
	    if( exists $rest{$key} ) {
		$name = $rest{$key};
		last;
	    }
	}
	if( ! defined $name ) { 
	    die("cannot parse this GFF3\n");
	}
	$names{$rest{'ID'}} = $name;
	print join("\t", @row, sprintf("%s %s",$group_name,$name)),"\n";
    }
}
