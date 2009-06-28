#!/usr/bin/perl -w
use strict;

my $GroupTag  = 'GenePrediction';
# for converting GFF3 into GFF2 for use by GLEAN
my %lookup;
while (<>) {
  next if /^\#/;
  chomp;
  my @row = split(/\t/,$_);
  my %lastcol= map { split(/=/, $_) } split(/;/,pop @row);
  if ( $row[2] eq 'mRNA') {
    my ($id,$name);
    for my $idtag ( qw(ID Parent Transcript) ) {
      if ( exists $lastcol{$idtag} ) {
	$id = $lastcol{$idtag};
	last;
      }
    }
    if ( ! defined $id ) {
      die("Unable to parse a valid ID or Parent field from these tags ", 
	  join(",", keys %lastcol), "\n");
    }
    if ( exists $lastcol{'Name'} ) {
      $name = $lastcol{'Name'};
    } else {
      $name = $id;
    }
    $lookup{$id} = $name;
    print join("\t", @row, sprintf("%s %s", $GroupTag, $name)),"\n";
  } elsif ( $row[2] eq 'CDS' ) {
    my ($id,$name);
    for my $idtag ( qw(Parent Transcript) ) {
      if ( exists $lastcol{$idtag} ) {
	$id = $lastcol{$idtag};
	last;
      }
    }
    if ( ! defined $id ) {
      die("No valid Parent field from these tags ",
	  join(",", keys %lastcol), "\n");
    }
    if ( exists $lookup{$id} ) {
      $name = $lookup{$id};
    } else {
      warn("cannot find parent name, using tag value as is\n");
      $name = $id;
    }
    print join("\t", @row, sprintf("%s %s", $GroupTag, $name)),"\n";
  }
}
