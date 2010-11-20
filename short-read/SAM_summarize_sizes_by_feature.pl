#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
use List::Util qw(sum);
my %expected_bases = map { $_ => 1 } qw(C A G T);
my @bases   = sort keys %expected_bases;
# provide SAM file(s) and feature file to process
my $grouping= 'type';

my ($genome,$features);
GetOptions(
	   'g|grouping:s'    => \$grouping,
	   'f|gff|feature:s' => \$features,
	   'genome|dna|fa:s' => \$genome,
	   );
die("must have feature file") unless defined $features && -f $features;
die("must have genome fasta") unless defined $genome && -f $genome;

my @bams = map { [$_ =~ /(\S+)\.bam/, Bio::DB::Sam->new(-bam => $_, -fasta => $genome)] } @ARGV;

open(my $fh => $features ) || die "$features: $!";

my %collected;
my $i = 0;
my %g;
while(<$fh>) {
    next if /^\#/;
    #next if /Simple_repeat|Low_complexity/;
    chomp;
    my @row = split(/\t/,$_);
    my $grp;
    my %group = map { split(/=/,$_) } split(';',pop @row);
    if( $grouping =~ /Note/i ) {	
	$grp =  $group{'Note'} || $group{'note'};
    } elsif( $grouping =~ /Name/i) {
    	$grp =  $group{'Name'} || $group{'name'};
    } elsif( $grouping =~ /ID/i ) {
	$grp =  $group{'ID'} || $group{'Id'};
    } elsif( $grouping =~ /type/i ) {
	# use the feature type as the grouping -- 3rd column
	$grp = $row[2];
    }
    $grp =~ s/\//_/g;
    $g{$grp}++;

    for my $bamd ( @bams ) {
	my ($name,$db) = @$bamd;
	my @alignments = $db->get_features_by_location(-seq_id => $row[0],
						       -start  => $row[3], 
						       -end    => $row[4]);
	for my $aln ( @alignments ) {
	    my $dna = $aln->query->dna ;
	    my %f;
	    for ( split('',$dna) ) {
		$f{$_}++;
	    }
	    next if keys %f <= 2; # drop those AAA or TTT runs
	    my $five_base = uc substr($dna,0,1);
	    next unless $expected_bases{$five_base};
	    my $len = $aln->length;
	    $collected{$name}->{$grp}->{$len}->{5}->{$five_base}++;	    
	    

	    my $three_base = uc substr($dna,-1,1);
	    next unless $expected_bases{$three_base};
	    $collected{$name}->{$grp}->{$len}->{3}->{$three_base}++;	    
	    
	}
    }
}
for my $r ( sort { $g{$b} <=> $g{$a} } keys %g ) {
    print join("\t",$r, $g{$r}),"\n";
}

my $odir = 'SAM_sizes_by_feature';
mkdir($odir);
open(my $R => ">$odir/summary_byfeature_barplot.R" ) || die $!;
while( my ($dbname,$coll) = each %collected ) {
    while( my ($type,$cts) = each %$coll ) {
	my %counts = %$cts;
	open(my $rpt => ">$odir/$dbname.$type.5_summary_sizes") || die $!;
	open(my $rptpct => ">$odir/$dbname.$type.5_summary_sizes.percent") || die $!;
	print $R "size5 <- read.table(\"$dbname.$type.5_summary_sizes\",header=T,sep=\"\\t\",row.names=1)\n";
	print $R "barplot(t(size5),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
    print $R "size5p <- read.table(\"$dbname.$type.5_summary_sizes.percent\",header=T,sep=\"\\t\",row.names=1)\n";
    print $R "barplot(t(size5p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname Reads mapped by size (percent) - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5p\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";

	my @lengths = sort { $a <=> $b } keys %counts;
	print $rpt join("\t", qw(LENGTH),  @bases),"\n";
	print $rptpct join("\t", qw(LENGTH), @bases),"\n";
	for my $c ( @lengths ) {
	    print $rpt join("\t", $c, map { $counts{$c}->{5}->{$_} || 0} @bases), "\n";
	    my $sum = sum ( map { $counts{$c}->{5}->{$_} || 0} @bases );
	    print $rptpct join("\t", $c, map { sprintf("%.2f",100*($counts{$c}->{5}->{$_} || 0)/$sum) } @bases),
	    "\n";
	}
	close($rpt);
	close($rptpct);
	if( 0 ) {
	    open($rpt => ">$odir/$dbname.$type.3_summary_sizes" ) || die $!;
	    open($rptpct => ">$odir/$dbname.$type..3_summary_sizes.percent" ) || die $!;
	    print $R "size3 <- read.table(\"$dbname.$type.3_summary_sizes\",header=T,sep=\"T\",row.names=1)\n";
	    print $R "barplot(t(size3),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size - 3' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size3\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
		print $R "size3p <- read.table(\"$dbname.$type.3_summary_sizes.percent\",header=T,sep=\"T\",row.names=1)\n";
	    print $R "barplot(t(size3p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size (percent) - 3' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size3p\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
		
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
    }
}
