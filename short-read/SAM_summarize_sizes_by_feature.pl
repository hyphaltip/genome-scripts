#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
use List::Util qw(sum);
my %expected_bases = map { $_ => 1 } qw(C A G T);
my @bases   = sort keys %expected_bases;

# provide SAM file(s) and feature file to process
my $grouping= 'type'; # which column to use for the feature type identification
my $minsize = 18;
my $maxsize = 30;
my $odir = 'SAM_sizes_by_feature';
my ($genome,$features,$debug);
GetOptions(
	   'v|verbose|debug!' => \$debug,
	   'min:i'            => \$minsize,
	   'max:i'            => \$maxsize,
	   'grouping:s'       => \$grouping,
	   'f|gff|feature:s'  => \$features,
	   'g|genome|dna|fa:s'=> \$genome,
	   'd|dir:s'          => \$odir
	   );
die("must have feature file") unless defined $features && -f $features;
die("must have genome fasta") unless defined $genome && -f $genome;
my @bams = map { [$_ =~ /^([^\.]+)\./, Bio::DB::Sam->new(-bam => $_, -fasta => $genome)] } @ARGV;

open(my $fh => $features ) || die "$features: $!";

my %collected;
my %by_chrom;
my $i = 0;
my %g;
while(<$fh>) {
    next if /^\#/;   
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
    } else {
	die("no defined field for group\n");
    }
    $grp =~ s/\//_/g;
    $g{$grp}++;
    for my $bamd ( @bams ) {
	my ($name,$db) = @$bamd;
	next if $debug && $db->length($row[0]) > 500000;
	my $segment = $db->segment($row[0], $row[3] => $row[4]);

	my $iterator = $segment->features(-iterator => 1);

	while (my $aln = $iterator->next_seq) {			
	    my $dna = $aln->query->dna ;
	    my %f;
	    for ( split('',$dna) ) { $f{$_}++ }
	    next if keys %f <= 2; # drop those AAA or TTT runs
	    
	    my $five_base = uc substr($dna,0,1);
	    next unless $expected_bases{$five_base};
	    my $len = $aln->length;
	    next if $len < $minsize || $len > $maxsize;
		
	    $collected{$name}->{$grp}->{$len}->{all}++;	    
	    $by_chrom{$name}->{$row[0]}->{$len}->{$grp}++;

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


mkdir($odir);
open(my $R => ">$odir/summary_byfeature_barplot.R" ) || die $!;

while( my ($dbname,$coll) = each %collected ) {
    print $R "pdf(\"$dbname\_reads_by_chrom_features_and_size.pdf\")\n";
    my @kinds = sort keys %$coll;    
    {
	open(my $rpt => ">$odir/$dbname.size_features_by_chrom.dat") || die $!;
	open(my $rptpct => ">$odir/$dbname.size_features_by_chrom_percent.dat") || die $!;
	
	print $rpt join("\t", qw(CHROM LENGTH), @kinds),"\n";
	print $rptpct join("\t", qw(CHROM LENGTH),@kinds),"\n";
	print $R "allsizechrom <- read.table(\"$dbname.size_features_by_chrom.dat\",header=T,sep=\"\\t\",row.names=NULL)\n";
	print $R "allsizechromp <- read.table(\"$dbname.size_features_by_chrom_percent.dat\",header=T,sep=\"\\t\",row.names=NULL)\n";
	for my $chrom ( sort keys %{$by_chrom{$dbname}} ) {
	    #open(my $rpt => ">$odir/$dbname.size_features_by_chrom.$chrom.cdat") || die $!;
	    #open(my $rptpct => ">$odir/$dbname.size_features_by_chrom_percent.$chrom.cdat") || die $!;
	
#	    print $rpt join("\t", qw(LENGTH), @kinds),"\n";
#	    print $rptpct join("\t", qw(LENGTH),@kinds),"\n";
	    #print $R "chrom <- read.table(\"$dbname.size_features_by_chrom.$chrom.cdat\",header=T,sep=\"\\t\",row.names=1)\n";
	    #print $R "chromp <- read.table(\"$dbname.size_features_by_chrom_percent.$chrom.cdat\",header=T,sep=\"\\t\",row.names=1)\n";	    	   

	    print $R "chrom <- subset(allsizechrom,allsizechrom\$CHROM == \"$chrom\");\n";
	    print $R "rownames(chrom) <- chrom[,2];\n";
	    print $R "chrom <- chrom[,3:",scalar @kinds+2,"];\n";

	    print $R "chromp <- subset(allsizechromp,allsizechromp\$CHROM == \"$chrom\");\n";
	    print $R "rownames(chromp) <- chromp[,2];\n";
	    print $R "chromp <- chromp[,3:",scalar @kinds+2,"];\n";
	    
	    print $R "barplot(t(chrom),xlab=\"Read Length\", ylab=\"Total \# Reads\", main=\"$chrom $dbname - Size and feature\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=chrom\$V1,legend=T,col=rainbow(",scalar @kinds+1,",start=0.1, end=.91),beside=F)\n";
		
	    print $R "barplot(t(chromp),xlab=\"Read Length\", ylab=\"Total \# Reads\", main=\" $chrom $dbname - Size and feature (percent)\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=chromp\$V1,legend=T,col=rainbow(",scalar @kinds,",start=0.1, end=.91),beside=F)\n";
	    for my $size ( sort { $a <=> $b } keys %{$by_chrom{$dbname}->{$chrom}} ) {
#		print $rpt join("\t", $size, map { $by_chrom{$dbname}->{$chrom}->{$size}->{$_} || 0 } @kinds ), "\n";
		print $rpt join("\t", $chrom, $size, map { $by_chrom{$dbname}->{$chrom}->{$size}->{$_} || 0 } @kinds ), "\n";
		my $sum = sum ( map { $by_chrom{$dbname}->{$chrom}->{$size}->{$_} || 0 } @kinds);
#		print $rptpct join("\t", $size, map { sprintf("%.2f",100 * ($by_chrom{$dbname}->{$chrom}->{$size}->{$_} || 0)/$sum) } @kinds), "\n";

		print $rptpct join("\t", $chrom, $size, 
				   map { sprintf("%.2f",100 * ($by_chrom{$dbname}->{$chrom}->{$size}->{$_} || 0)/$sum) } @kinds), "\n";
	    }
	}
    }

    my %size2ftype_count;
    print $R "pdf(\"$dbname\_reads_by_feature.pdf\");\n";
    while( my ($type,$cts) = each %$coll ) {
	my %counts = %$cts;
	open(my $rpt => ">$odir/$dbname.$type.5_summary_sizes") || die $!;
	open(my $rptpct => ">$odir/$dbname.$type.5_summary_sizes.percent") || die $!;
	print $R "size5 <- read.table(\"$dbname.$type.5_summary_sizes\",header=T,sep=\"\\t\",row.names=1)\n";
	print $R "barplot(t(size5),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5\$V1,legend=T,col=rainbow(5,start=0.1, end=.91),beside=F)\n";

	    print $R "size5p <- read.table(\"$dbname.$type.5_summary_sizes.percent\",header=T,sep=\"\\t\",row.names=1)\n";
	print $R "barplot(t(size5p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size (percent) - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5p\$V1,legend=T,col=rainbow(5,start=0.1, end=.91),beside=F)\n";

	    my @lengths = sort { $a <=> $b } keys %counts;
	print $rpt join("\t", qw(LENGTH),  @bases),"\n";
	print $rptpct join("\t", qw(LENGTH), @bases),"\n";
	for my $c ( @lengths ) {
	    print $rpt join("\t", $c, map { $counts{$c}->{5}->{$_} || 0} @bases), "\n";
	    my $sum = sum ( map { $counts{$c}->{5}->{$_} || 0} @bases );
	    $size2ftype_count{$c}->{$type} = $sum;

	    print $rptpct join("\t", $c, map { sprintf("%.2f",100*($counts{$c}->{5}->{$_} || 0)/$sum) } @bases),
	    "\n";
	}
	close($rpt);
	close($rptpct);    	
	
#	if( 0 ) {
#	    open($rpt => ">$odir/$dbname.$type.3_summary_sizes" ) || die $!;
#	    open($rptpct => ">$odir/$dbname.$type..3_summary_sizes.percent" ) || die $!;
#	    print $R "size3 <- read.table(\"$dbname.$type.3_summary_sizes\",header=T,sep=\"T\",row.names=1)\n";
#	    print $R "barplot(t(size3),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size - 3' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size3\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
#		print $R "size3p <- read.table(\"$dbname.$type.3_summary_sizes.percent\",header=T,sep=\"T\",row.names=1)\n";
#	    print $R "barplot(t(size3p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size (percent) - 3' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size3p\$V1,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
#		
#		@lengths = sort { $a <=> $b } keys %counts;
#	    @bases   = sort keys %expected_bases;
#	    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
#	    print $rptpct join("\t", qw(LENGTH), @bases),"\n";
#
#	    for my $c ( @lengths ) {
#		print $rpt join("\t", $c, map { $counts{$c}->{3}->{$_} } @bases), "\n";
#		my $sum = sum ( map { $counts{$c}->{3}->{$_} } @bases );
#		print $rptpct join("\t", $c, map { sprintf("%.2f",100*$counts{$c}->{3}->{$_}/$sum) } @bases),
#		"\n";
#	    }		
    }
    open(my $rpt => ">$odir/$dbname.feature_sizes.dat") || die $!;
    open(my $rptpct => ">$odir/$dbname.feature_sizes_percent.dat") || die $!;    
    print $rpt join("\t", qw(LENGTH), @kinds), "\n";
    print $rptpct join("\t", qw(LENGTH), @kinds), "\n";
    for my $size ( sort { $a <=> $b } keys %size2ftype_count ) {
	print $rpt join("\t", $size, map { $size2ftype_count{$size}->{$_} || 0 } @kinds),"\n";
	my $sum = sum ( map { $size2ftype_count{$size}->{$_} || 0 } @kinds );
	print $rptpct join("\t", $size, map { sprintf("%.2f",100 * ($size2ftype_count{$size}->{$_} || 0) / $sum) } @kinds),"\n";
    }
    close($rpt);
    close($rptpct);
    print $R "pdf(\"$dbname\_reads_by_feature_and_size_combined.pdf\")\n";
    print $R "sizef <- read.table(\"$dbname.feature_sizes.dat\",header=T,sep=\"\\t\",row.names=1)\n";
    print $R "sizefp <- read.table(\"$dbname.feature_sizes_percent.dat\",header=T,sep=\"\\t\",row.names=1)\n";

    print $R "barplot(t(sizef),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname Reads mapped by size and feature type\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=sizef\$V1,legend=T,col=rainbow(",scalar @kinds,",start=0.1, end=.91),beside=F)\n";
    print $R "barplot(t(sizefp),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname Reads mapped by size (percent) and feature type\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=sizefp\$V1,legend=T,col=rainbow(",scalar @kinds,",start=0.1, end=.91),beside=F)\n";
    
}
