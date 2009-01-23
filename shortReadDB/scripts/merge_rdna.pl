#!/usr/bin/perl -w
use strict;
use List::Util qw(min max);

# seqid, source, type, start, end, score, strand, frame, group
#  0       1       2    3      4    5      6       7,      8

print "##gff-version 3\n";

my %groups;
my %exons;
my @order;
while(<>) {
    next if /^\#/;
    chomp;
    my $line = $_;
    my @line = split(/\t/,$_);
    my %col = map { split(/=/,$_) } pop @line;
    my ($id);
    for my $t ( qw(ID Parent Transcript) ) {
	next unless exists $col{$t};
	$id = $col{$t};
	unless( exists $groups{$id} ) {
	    push @order, $id;
	}
	last;
    }
    if( ! defined $id ) {
	warn("cannot get an ID from $line\n");
	next;
    }
    $groups{$id}->[0] =  $line[0];
    $line[1] = 'BLAST_rDNA';
    $groups{$id}->[1] =  $line[1];
    $groups{$id}->[2] =  'rDNA';

    $groups{$id}->[3] = defined $groups{$id}->[3] ? min($groups{$id}->[3],$line[3]): $line[3];
    $groups{$id}->[4] = defined $groups{$id}->[4] ? max($groups{$id}->[4],$line[4]): $line[4];
    $groups{$id}->[5] =  $line[5];
    $groups{$id}->[6] =  $line[6];
    $groups{$id}->[7] =  $line[7];
    $groups{$id}->[8] = sprintf("ID=%s;Name=%s",$id,$id);
#    push @{$exons{$id}}, join("\t",@line, 
#			      sprintf("ID=%s.e%d;Parent=%s;Name=%s.e%d",$id,scalar @{$exons{$id} || []}, $id,
#				      $id, scalar @{$exons{$id} || []}));
}

for my $f ( @order ) {
    print join("\t",@{$groups{$f}}),"\n";
#    print join("\n",@{$exons{$f}}),"\n";
}
