#!/usr/bin/perl -w
use strict;

# this script expects a pfam folder to contain either
# hmmscan output with --domtblout 
# OR 
# hmmpfam output converted to tab deliminted with hmmer_to_table.pl
# The files are expected to end in .tab though that is modifiable
# with -ext cmdline option
# If you use HMMER2 output you need to specify this with
# -hmmer 2 on the cmdline (-hmmer 3 is the default)
#
# the OrthoMCL 2 output file from orthoMclToGroups output is
# provided as the remaining command line argument
# 
# by default domains with evalue > 0.01 are discarded if you want to adjust
# this cutoff use -e or --evalue option

use Getopt::Long;
my $pfamdir;
my $pfamext = 'tab';
my $hmmerversion = 3;
my $evalue_cutoff = 0.01;
GetOptions(
	   'd|p|pfam:s'  => \$pfamdir,
	   'ext:s'     => \$pfamext,
	   'hmmer|hv:s'=> \$hmmerversion,
	   'e|evalue:s' => \$evalue_cutoff,
	   );

die "must provide a folder with pfam table data (either domtblout or hmmer2table output\n" unless $pfamdir && -d $pfamdir;

my %pfam;
opendir(PFAM, $pfamdir) || die "cannot open $pfamdir: $!";
for my $file ( readdir(PFAM) ) {
    next unless ( $file =~ /\.\Q$pfamext\E$/);
    open(my $fh => "$pfamdir/$file" ) || die $!;
    if( $hmmerversion == 3 ) { #parse HMMER3 domtblout output
	while(<$fh>) {
	    next if /^\#/;
	    my ($domain,$acesssion,$tlen, $qname, $qacc,
		$qlen, $seq_evalue, $seq_score,$seq_bias,
		$n, $totalhits, $dom_cvalue, $dom_ievalue,$dom_score,
		$dom_bias) = split(/\s+/,$_);
	    $pfam{$qname}->{$domain}++ if $dom_cvalue < $evalue_cutoff;
	    $qname =~ s/:/\|/;
	    $pfam{$qname}->{$domain}++ if $dom_cvalue < $evalue_cutoff;
	}
    } elsif( $hmmerversion ==2 ) { # parse HMMER2 (hmmer_to_table output
       while(<$fh>) {
	   next if /^\#/;
	   my ($qname, $qstart, $qend, $domain,
	       $domstart, $domend, $score, $evalue) = split;	       
	   $pfam{$qname}->{$domain}++ if $evalue < $evalue_cutoff;
       }
   } else { 
       warn("unknown HMMER version (2 or 3 are expected)\n");
   }
}

my %sp;
my @data;
while(<>) {
    my ($group,@orthologs) = split;
    $group =~ s/://;
    my %count;
    my %domains;	
    for my $gene ( @orthologs ) {
	my ($sp)= split(/\|/,$gene);
	$count{$sp}++;
	$sp{$sp} = 1; # collecting all the species names
	while( my ($dom,$count) = each %{$pfam{$gene} || {}} ) {
	    $domains{$dom} += $count;
	}
    }    
    push @data,[ $group, \%count, \%domains];
}
my @sp = sort keys %sp;
print join("\t", 'FAMILY', @sp, 'DOMAINS'),"\n";
for my $d ( @data ) {
    my ($group,$count,$domains) = @$d;
    print join("\t", $group, 
	       (map { $count->{$_} || 0 } @sp),
	       join(",", map { sprintf("%s (%d)",$_, 
				       $domains->{$_}) } 
		    sort { $domains->{$b} <=> $domains->{$a} } 
		    keys %$domains)),"\n";    
}
