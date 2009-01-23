#!/usr/bin/perl -w
use Getopt::Long;
use strict;
my %count;
my $i;
my $debug;
GetOptions('v|verbose|debug!' => \$debug);
while(<>) {
	if( /^\#/ ) {
		print;
		next;
	}
	chomp;
	my @line = split(/\t/,$_);
	if( $line[2] =~ /^(CDS|exon)/ ) {
		$line[-1] =~ s/\.gene//; # fix
	}
	warn("last is $line[-1]\n") if $debug;
	my %group = map { split(/=/,$_) } split (/;/,pop @line);
	unless( $group{ID} ) {	
	my $stem = sprintf("%s.%s",$group{Parent},$line[2]);
		$group{ID} = sprintf("%s.%d",$stem, ++$count{$stem});
	}
	$group{Name} = $group{ID};
	push @line, join(";",map { sprintf("%s=%s",$_, $group{$_}) } grep { exists $group{$_} } qw(ID Parent Name Alias));
 	print join("\t", @line),"\n";	
	last if $i++ > 10 && $debug;
}
