#!/usr/bin/perl
use strict;
use warnings;
my %count;
while(<>) {
    next if /^\s+$/;
    if( /^#/ ) {
	print $_;
	next;
    }
    chomp;
    my @row = split(/\t/,$_);
    my @order;
    my %grp = map { 	
	my ($key,$val) = split (/=/, $_);
	push @order, $key;
	($key,$val); } split(/;/,pop @row);

    if( $row[2] eq 'CDS' ) {
	$grp{'ID'} .= ".cds".++$count{$grp{'Parent'}};
    }
    my $ninth = join(";",map { sprintf("%s=%s",$_,$grp{$_}) } @order);
    print join("\t", @row,$ninth),"\n";
}
