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
my $features;
my $debug = 0;
GetOptions('v|verbose|debug!' => \$debug,
	   'd|dir:s'   => \$dir,
	   'min:i'     => \$minsize,
	   'max:i'     => \$maxsize,
	   'f|skip:s'  => \$features,
	   'g|genome:s'=> \$genome);

mkdir($dir) unless -d $dir;

my %skip;
if( $features ) {

    open(my $fh => $features ) || die "$features: $!";
    while(<$fh>) {
	next if /^\#/;
	chomp;
	my ($seqid,$src,$type,$start,$end,undef,$strand,@rest) = split;

	push @{$skip{$seqid}}, [ (sort { $a <=> $b} ( $start, $end)) ,
				 $strand eq '+' ? 1 : -1, pop @rest];
    }
}
open(my $R => ">$dir/summary_shortread_barplot.R" ) || die $!;
for my $file ( @ARGV ) {
    my $base = $file;
    $base =~ s/\.bam$//;
    print $R "pdf(\"$base.pdf\")\n";
    $base =~ s/\.bam$//;
    my $db = Bio::DB::Sam->new(-bam => $file,
			       -fasta=>$genome);
    my %counts;
    for my $id ( $db->seq_ids ) {
	my $start_seq = 1;
	my $length = $db->length($id);
	for my $skipfeature ( sort { $a->[0] <=> $b->[0] } @{$skip{$id} || []} ) {
	    my ($sfeat_start,$sfeat_end, $sfeat_strand,$name) = @$skipfeature;
	    if( ($sfeat_start-1) <= $start_seq ) {
		$start_seq = $sfeat_end + 1;		
		next;
	    }
	    my $segment = $db->segment($id,$start_seq => $sfeat_start - 1);
	    $start_seq = $sfeat_end + 1;
	    my $iterator     = $segment->features(-iterator=>1);
	    
	    while (my $align = $iterator->next_seq) { 			    
		my $seq = $align->query->dna;
		
		my %f;
		for ( split('',$seq) ) { $f{$_}++ }
		next if keys %f <= 2; # drop those AAA or TTT runs
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
	if( $start_seq < $length ) {
	    my $segment = $db->segment($id,$start_seq => $length);
	    my $iterator     = $segment->features(-iterator=>1);
	    
	    while (my $align = $iterator->next_seq) { 			    
		my $seq = $align->query->dna;		
		my %f;
		for ( split('',$seq) ) { $f{$_}++ }
		next if keys %f <= 2; # drop those AAA or TTT runs
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
	last if $debug;
    }
    
    open(my $rpt => ">$dir/$base.5_summary_sizes" ) || die $!;
    open(my $rptpct => ">$dir/$base.5_summary_sizes.percent" ) || die $!;
    print $R "size5 <- read.table(\"$base.5_summary_sizes\",header=T,sep=\"\\t\",row.names=1)\n";
    print $R "barplot(t(size5),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size - 5' base\",",
    "legend=T,col=c(\"red\",\"blue\",\"green\",\"orange\"),beside=T)\n";
    print $R "size5p <- read.table(\"$base.5_summary_sizes.percent\",header=T,sep=\"\\t\",row.names=1)\n";

# add in total percentage plots
    print $R "barplot(t(size5p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size (percent) - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";

    my @lengths = sort { $a <=> $b } keys %counts;    		

    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
    print $rptpct join("\t", qw(LENGTH), @bases),"\n";
    for my $c ( @lengths ) {
	print $rpt join("\t", $c, map { $counts{$c}->{5}->{$_} || 0 } @bases), "\n";
	my $sum = sum ( map { $counts{$c}->{5}->{$_} || 0 } @bases );
	print $rptpct join("\t", $c, map { sprintf("%.2f",100*($counts{$c}->{5}->{$_} || 0)/$sum) } @bases),
	"\n";
    }
    close($rpt);
    close($rptpct);

    open($rpt => ">$dir/$base.3_summary_sizes" ) || die $!;
    open($rptpct => ">$dir/$base.3_summary_sizes.percent" ) || die $!;
    print $R "size3 <- read.table(\"$base.3_summary_sizes\",header=T,sep=\"\\t\",row.names=1)\n";
    print $R "barplot(t(size3),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size - 3' base\",",
    "legend=T,col=c(\"red\",\"blue\",\"green\",\"orange\"),beside=T)\n";
    print $R "size3p <- read.table(\"$base.3_summary_sizes.percent\",header=T,sep=\"\\t\",row.names=1)\n";

# add in total percentage plots
    print $R "barplot(t(size3p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size (percent) - 3' base\",",
    "legend=T,col=c(\"red\",\"blue\",\"green\",\"orange\"),beside=T)\n";
    
    @lengths = sort { $a <=> $b } keys %counts;
    @bases   = sort keys %expected_bases;
    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
    print $rptpct join("\t", qw(LENGTH), @bases),"\n";

    for my $c ( @lengths ) {
	print $rpt join("\t", $c, map { $counts{$c}->{3}->{$_} || 0} @bases), "\n";
	my $sum = sum ( map { $counts{$c}->{3}->{$_} || 0} @bases );
	print $rptpct join("\t", $c, map { sprintf("%.2f",100*($counts{$c}->{3}->{$_} || 0)/$sum) } @bases),
	"\n";
    }

}
