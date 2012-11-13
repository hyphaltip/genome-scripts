#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use List::Util qw(sum);
use Bio::DB::Sam;

# this expects BAM files as the ARGV
my %expected_bases = map { $_ => 1 } qw(C A G T);
my @bases   = sort keys %expected_bases;
my $dir = 'size_summary_stranded';
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
	    $start_seq  = $sfeat_end + 1;
	    my $iterator = $segment->features(-iterator=>1);
	    
	    while (my $align = $iterator->next_seq) { 			    
		my $seq = $align->query->dna;
		my $strand = $align->strand < 0 ? 'minus' : 'plus';
		my %f;
		for ( split('',$seq) ) { $f{$_}++ }
		next if keys %f <= 2; # drop those AAA or TTT runs
		my $length = $align->query->length;
		next if $length < $minsize || $length > $maxsize;
		my $five_base = uc substr($seq,0,1); # get 5' base
		my $three_base = uc substr($seq,-1,1); # get 3' base
		next unless $expected_bases{$five_base};
		$counts{$length}->{$strand}->{5}->{$five_base}++;
		next unless $expected_bases{$three_base};
		$counts{$length}->{$strand}->{3}->{$three_base}++;
	    }
	}
	if( $start_seq < $length ) {
	    my $segment = $db->segment($id,$start_seq => $length);
	    my $iterator     = $segment->features(-iterator=>1);
	    
	    while (my $align = $iterator->next_seq) { 			    
		my $seq = $align->query->dna;		
		my $strand = $align->strand < 0 ? 'minus' : 'plus';
		my %f;
		for ( split('',$seq) ) { $f{$_}++ }
		next if keys %f <= 2; # drop those AAA or TTT runs
		my $length = $align->query->length;
		next if $length < $minsize || $length > $maxsize;
		my $five_base = uc substr($seq,0,1); # get 5' base
		my $three_base = uc substr($seq,-1,1); # get 3' base
		next unless $expected_bases{$five_base};
#		$counts{$length}->{5}->{$five_base}++;
		$counts{$length}->{$strand}->{5}->{$five_base}++;
		next unless $expected_bases{$three_base};
#		$counts{$length}->{3}->{$three_base}++;
		$counts{$length}->{$strand}->{3}->{$three_base}++;
	    }
	}
	last if $debug;
    }
    for my $strand ( qw(plus minus) ) {
	for my $readend ( 5,3) {
	    open(my $rpt => ">$dir/$base.$readend\_summary_sizes_plus" ) || die $!;
	    open(my $rptpct => ">$dir/$base.$readend\_summary_sizes.percent_plus" ) || die $!;
	    print $R "size <- read.table(\"$base.$readend\_summary_sizes_plus\",header=T,sep=\"\\t\",row.names=1)\n";
	    print $R "barplot(t(size),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size - $readend' base $strand strand\",",
	    "legend=T,col=c(\"red\",\"blue\",\"green\",\"orange\"),beside=T)\n";
	    print $R "sizep <- read.table(\"$base.$readend\_summary_sizes.percent_plus\",header=T,sep=\"\\t\",row.names=1)\n";
	    
# add in total percentage plots
	    print $R "barplot(t(sizep),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"Reads mapped by size (percent) - $readend' base $strand strand\",space=0.1,cex.axis=0.8,las=1,cex=0.8,legend=T,col=rainbow(5,start=.1, end=.91),beside=F)\n";
	    
	    my @lengths = sort { $a <=> $b } keys %counts;    		
	    
    print $rpt join("\t", qw(LENGTH),  @bases),"\n";
	    print $rptpct join("\t", qw(LENGTH), @bases),"\n";
	    for my $c ( @lengths ) {
		print $rpt join("\t", $c, map { $counts{$c}->{$strand}->{$readend}->{$_} || 0 } @bases), "\n";
		my $sum = sum ( map { $counts{$c}->{$strand}->{$readend}->{$_} || 0 } @bases );
		print $rptpct join("\t", $c, map { sprintf("%.2f",100*($counts{$c}->{$strand}->{$readend}->{$_} || 0)/$sum) } @bases),
		"\n";
	    }
	    close($rpt);
	    close($rptpct);
	}
    } 
}
