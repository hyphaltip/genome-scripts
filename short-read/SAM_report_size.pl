#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use List::Util qw(sum);
use Bio::DB::Sam;

# this expects BAM files as the ARGV
my %expected_bases = map { $_ => 1 } qw(C A G T);
my @bases   = sort keys %expected_bases;
my $dir = 'size_summary';
my $minsize = 18;
my $maxsize = 30;
my $genome;
GetOptions('d|dir:s' => \$dir,
	   'min:i'     => \$minsize,
	   'max:i'     => \$maxsize,
	   'g|genome:s'=> \$genome);

mkdir($dir) unless -d $dir;

open(my $R => ">$dir/summary_shortread_barplot.R" ) || die $!;
for my $file ( @ARGV ) {
    my $base = $file;
    print $R "pdf(\"$base.pdf\")\n";
    $base =~ s/\.bam$//;
    my $db = Bio::DB::Sam->new(-bam => $file,
			       -fasta=>$genome);
    my %counts;
    for my $id ( $db->seq_ids ) {
	my $segment = $db->segment($id);
	my $iterator     = $segment->features(-iterator=>1);
	while (my $align = $iterator->next_seq) { 
	    my $seq = $align->query->dna;
	    my $length = $align->query->length;
	    next if $length < $minsize || $length > $maxsize;
	    my $five_base = uc substr($seq,0,1); # get 5' base
	    my $three_base = uc substr($seq,-1,1); # get 3' base
	    next unless $expected_bases{$five_base};
	    $counts{$length}->{5}->{$five_base}++;
	    next unless $expected_bases{$three_base};
	    $counts{$length}->{3}->{$three_base}++;
	}
    }
    
    open(my $rpt => ">$dir/$base.5_summary_sizes" ) || die $!;
    open(my $rptpct => ">$dir/$base.5_summary_sizes.percent" ) || die $!;
    print $R "size5 <- read.table(\"$base.5_summary_sizes\",header=T,sep=\"\\t\",row.names=1)\n";
    print $R "barplot(t(size5),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
    print $R "size5p <- read.table(\"$base.5_summary_sizes.percent\",header=T,sep=\"\\t\",row.names=1)\n";
    print $R "barplot(t(size5p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size (percent) - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5p\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";

    my @lengths = sort { $a <=> $b } keys %counts;    
    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
    print $rptpct join("\t", qw(LENGTH), @bases),"\n";
    for my $c ( @lengths ) {
	print $rpt join("\t", $c, map { $counts{$c}->{5}->{$_} } @bases), "\n";
	my $sum = sum ( map { $counts{$c}->{5}->{$_} } @bases );
	print $rptpct join("\t", $c, map { sprintf("%.2f",100*$counts{$c}->{5}->{$_}/$sum) } @bases),
	"\n";
    }
    close($rpt);
    close($rptpct);

    open($rpt => ">$dir/$base.3_summary_sizes" ) || die $!;
    open($rptpct => ">$dir/$base.3_summary_sizes.percent" ) || die $!;
    print $R "size3 <- read.table(\"$base.3_summary_sizes\",header=T,sep=\"T\",row.names=1)\n";
    print $R "barplot(t(size3),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size - 3' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size3\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
    print $R "size3p <- read.table(\"$base.3_summary_sizes.percent\",header=T,sep=\"T\",row.names=1)\n";
    print $R "barplot(t(size3p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size (percent) - 3' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size3p\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
    
    @lengths = sort { $a <=> $b } keys %counts;
    @bases   = sort keys %expected_bases;
    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
    print $rptpct join("\t", qw(LENGTH), @bases),"\n";

    for my $c ( @lengths ) {
	print $rpt join("\t", $c, map { $counts{$c}->{3}->{$_} } @bases), "\n";
	my $sum = sum ( map { $counts{$c}->{3}->{$_} } @bases );
	print $rptpct join("\t", $c, map { sprintf("%.2f",100*$counts{$c}->{3}->{$_}/$sum) } @bases),
	"\n";
    }

}

