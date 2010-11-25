#!/usr/bin/perl -w
use strict;

# This script is run on files produced by 
# SAM_summarize_fwindow.pl

my $min_total = 100;
my %combo;
for my $file ( @ARGV ) {
    next unless ( $file =~ /(\S+)\.counts\.dat$/);
    my $stem = $1;
    open(my $fh => $file)|| die $!;
    my $hdr = <$fh>;
    my $i = 0;
    chomp $hdr;
    my %hdr = map { $_ => $i++ } split(/\t/,$hdr);
    while(<$fh>) {
	chomp;
	my @row = split(/\t/,$_);
	my ($f,$length,$same,$opp) = 
	    map { $row[ $hdr{$_} ] } qw(FEATURE LENGTH SAME OPP);
	$combo{$f}->{$stem}->{$length} = $same + $opp;
    }    
}
open(my $R => ">plot_fcount.R") || die $!;
print $R "pdf(\"plot_fcount.pdf\");\n";
for my $feat ( sort keys %combo ) {
    open(my $fh => ">$feat.tab") || die $!;
    my %lens;
    my $total = 0;
    for my $stem ( sort keys %{$combo{$feat}} ) {
	for my $l ( keys %{$combo{$feat}->{$stem}} ) {
	    $lens{$l}++;
	    $total += $combo{$feat}->{$stem}->{$l};
	}
    }
    next if $total < $min_total;
    my @lenslst = sort { $a <=> $b} keys %lens;
    print $fh join("\t", 'LIBRARY',@lenslst), "\n";
    for my $stem ( sort keys %{$combo{$feat}} ) {
	print $fh join("\t", $stem, map { $combo{$feat}->{$stem}->{$_} || 0 } @lenslst),"\n";
    }
    
    print $R "f$feat <- read.table(\"$feat\.tab\",sep=\"\\t\",header=T,row.names=1);\n";
    print $R "barplot(t(f$feat),xlab=\"Library\", ylab=\"Total Reads\", main=\"$feat - read size\",space=0.1,cex.axis=0.8,las=1,cex=0.8,legend=T,col=rainbow(",scalar @lenslst+1,",start=0.1, end=.91),beside=F)\n";
    print $R "fp$feat <- prop.table(as.matrix(f$feat),margin=1)*100\n";
    print $R "barplot(t(fp$feat),xlab=\"Library\", ylab=\"Total Reads\", main=\"$feat - read size \%\",space=0.1,cex.axis=0.8,las=1,cex=0.8,legend=T,col=rainbow(",scalar @lenslst+1,",start=0.1, end=.91),beside=F)\n";
}
