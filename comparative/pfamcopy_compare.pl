#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);

my $cutoff = 1e-3;
my $dir = shift || 'pfam';
my $ext = '\.domtbl\.tab';
opendir(DIR,$dir) || die "$dir: $!";
my %domains;
my %spnames;
for my $file ( readdir(DIR) ) {
    if( $file =~ /(\S+)$ext/ ) {
	my $sp = $1;
	$spnames{$sp}++;
#	print "$1 $file\n";
	open(my $fh => "$dir/$file" ) || die $!;
	while(<$fh>) {
	    next if /^\#/;
	    chomp;
	    my ($domain,$acc,$tlen,$gene,$trans,$qlen,$evalue_full,$score_full,$bias_full,
		$n, $tot,$c_evalue,$i_evalue, $domain_score,$domain_bias,
		$hmm_start,$hmm_end, $ali_start, $ali_end, 
		$env_start, $env_end,
		$accuracy, 
		$description) = split(/\s+/,$_,24);
	    next if $i_evalue > $cutoff;
	    $domains{$domain}->{$sp}++;
	}
    }
}
my @sp = sort keys %spnames;
print join("\t",
	   qw(DOMAIN TOTAL),
	   @sp),"\n";

for my $dom ( sort { $b->[1] <=> $a->[1] } 
	      map { [$_,sum(values %{$domains{$_}}), ] } 
	      keys %domains ) {
    print join("\t", $dom->[0], $dom->[1],
	       map { $domains{$dom->[0]}->{$_} || 0 } @sp), "\n";    
}
