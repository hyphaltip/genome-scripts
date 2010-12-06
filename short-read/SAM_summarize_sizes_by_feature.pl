#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
use File::Spec;
use List::Util qw(sum);
my %expected_bases = map { $_ => 1 } qw(C A G T);
my @bases   = sort keys %expected_bases;

# provide SAM file(s) and feature file to process
my $grouping= 'type'; # which column to use for the feature type identification
my $minsize = 18;
my $maxsize = 30;
my $odir = 'SAM_sizes_by_feature';
my ($genome,$features,$debug);
my $compute_nofeature = 0;
GetOptions(
	   'v|verbose|debug!' => \$debug,
	   'min:i'            => \$minsize,
	   'max:i'            => \$maxsize,
	   'grouping:s'       => \$grouping,
	   'f|gff|feature:s'  => \$features,
	   'g|genome|dna|fa:s'=> \$genome,
	   'nf|nofeature!'    => \$compute_nofeature,
	   'd|dir:s'          => \$odir
	   );
die("must have feature file") unless defined $features && -f $features;
die("must have genome fasta") unless defined $genome && -f $genome;
open(my $fh => $features ) || die "$features: $!";

my @bams = map { [$_ =~ /^([^\.]+)\./, Bio::DB::Sam->new(-bam => $_, -fasta => $genome)] } @ARGV;
my %chrom_lens;
for my $chrom ( $bams[0]->[1]->seq_ids ) {
    $chrom_lens{$chrom} = $bams[0]->[1]->length($chrom);
}

my (undef,undef,$ffilename) = File::Spec->splitpath($features);
$ffilename =~ s/\.gff\S+$//;
open(my $groupsfh => ">seen_groups.$ffilename.out");
my %collected;
my %by_chrom;
my $i = 0;
my %g;
my %genome_features;
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
	my ($start,$end) = sort { $a <=> $b} ($row[3],$row[4]);
	my $segment = $db->segment($row[0], $start => $end);
	push @{$genome_features{$row[0]}}, [$start => $end];
	
	my $iterator = $segment->features(-iterator => 1);

	while (my $aln = $iterator->next_seq) {			
	    my $len = $aln->length;
	    next if $len < $minsize || $len > $maxsize;
		
	    my $dna = $aln->query->dna ;
	    my %f;
	    for ( split('',$dna) ) { $f{$_}++ }
	    next if keys %f <= 2; # drop those AAA or TTT runs
	    
	    my $five_base = uc substr($dna,0,1);
	    my $three_base = uc substr($dna,-1,1);

	    # stay super strict here
	    next unless $expected_bases{$three_base} && $expected_bases{$five_base};
	    $collected{$name}->{$grp}->{$len}->{5}->{$five_base}++;	    	    
	    $collected{$name}->{$grp}->{$len}->{3}->{$three_base}++;	    

	    
	    $collected{$name}->{$grp}->{$len}->{all}++;	    
	    $by_chrom{$name}->{$row[0]}->{$len}->{$grp}++;
	}
    }
}

if( $compute_nofeature ) {
    my $hole_set = &compute_holes(\%genome_features, \%chrom_lens);

    while( my ($hole_chrom,$holes) = each %$hole_set ) {
	for my $hole ( @{$holes} ) {
	    for my $bamd ( @bams ) {
		my ($name,$db) = @$bamd;
		next if $debug && $db->length($hole_chrom) > 500000;
		my $segment = $db->segment($hole_chrom, $hole->[0] => $hole->[1]);

		my $iterator = $segment->features(-iterator => 1);	    
		while (my $aln = $iterator->next_seq) {			
		    my $len = $aln->length;
		    next if $len < $minsize || $len > $maxsize;

		    my $dna = $aln->query->dna ;
		    my %f;
		    for ( split('',$dna) ) { $f{$_}++ }
		    next if keys %f <= 2; # drop those AAA or TTT runs


		    my $five_base = uc substr($dna,0,1);
		    my $three_base = uc substr($dna,-1,1);

		    # stay super strict here
		    next unless $expected_bases{$three_base} && $expected_bases{$five_base};

		    $collected{$name}->{'nofeature'}->{$len}->{5}->{$five_base}++;	    	    


		    $collected{$name}->{'nofeature'}->{$len}->{3}->{$three_base}++;	    

		    $collected{$name}->{'nofeature'}->{$len}->{all}++;		
		    $by_chrom{$name}->{$hole_chrom}->{$len}->{'nofeature'}++;

		}
	    }
	}
    }
}
for my $r ( sort { $g{$b} <=> $g{$a} } keys %g ) {
    print $groupsfh join("\t",$r, $g{$r}),"\n";
}


mkdir($odir);
open(my $R => ">$odir/summary_byfeature_barplot.R" ) || die $!;

while( my ($dbname,$coll) = each %collected ) {
    print $R "pdf(\"$dbname\_reads_by_chrom_features_and_size.pdf\")\n";
    my @kinds = sort keys %$coll;    
    {
	open(my $rpt => ">$odir/$dbname.size_features_by_chrom.dat") || die $!;
	
	print $rpt join("\t", qw(CHROM LENGTH), @kinds),"\n";
	print $R "allsizechrom <- read.table(\"$dbname.size_features_by_chrom.dat\",header=T,sep=\"\\t\",row.names=NULL)\n";
	for my $chrom ( sort keys %{$by_chrom{$dbname}} ) {
	    print $R "chrom <- subset(allsizechrom,allsizechrom\$CHROM == \"$chrom\");\n";
	    print $R "rownames(chrom) <- chrom[,2];\n";
	    print $R "chrom <- chrom[,3:",scalar @kinds+2,"];\n";
	    
	    print $R "chromp <- prop.table(as.matrix(chrom),margin=1)*100;\n";
	    
	    
	    print $R "barplot(t(chrom),xlab=\"Read Length\", ylab=\"Total \# Reads\", main=\"$chrom $dbname - Size and feature\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=chrom\$V1,legend=T,col=rainbow(",scalar @kinds+1,",start=0.1, end=.91),beside=F)\n";
		
	    print $R "barplot(t(chromp),xlab=\"Read Length\", ylab=\"Total \# Reads\", main=\" $chrom $dbname - Size and feature (percent)\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=chrom\$V1,legend=T,col=rainbow(",scalar @kinds,",start=0.1, end=.91),beside=F)\n";
	    for my $size ( sort { $a <=> $b } keys %{$by_chrom{$dbname}->{$chrom}} ) {
		print $rpt join("\t", $chrom, $size, map { $by_chrom{$dbname}->{$chrom}->{$size}->{$_} || 0 } @kinds ), "\n";
	    }
	}
    }

    my %size2ftype_count;
    print $R "pdf(\"$dbname\_reads_by_feature.pdf\");\n";
    while( my ($type,$cts) = each %$coll ) {
	my %counts = %$cts;
	open(my $rpt => ">$odir/$dbname.$type.5_summary_sizes") || die $!;
	print $R "size5 <- read.table(\"$dbname.$type.5_summary_sizes\",header=T,sep=\"\\t\",row.names=1)\n";
	print $R "barplot(t(size5),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5\$V1,legend=T,col=rainbow(5,start=0.1, end=.91),beside=F)\n";

	    print $R "size5p <- prop.table(as.matrix(size5),margin=1)*100;\n";
	
	print $R "barplot(t(size5p),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname $type Reads mapped by size (percent) - 5' base\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=size5\$V1,legend=T,col=rainbow(5,start=0.1, end=.91),beside=F)\n";

	    my @lengths = sort { $a <=> $b } keys %counts;
	print $rpt join("\t", qw(LENGTH),  @bases),"\n";
	for my $c ( @lengths ) {
	    print $rpt join("\t", $c, map { $counts{$c}->{5}->{$_} || 0} @bases), "\n";
	    my $sum = sum ( map { $counts{$c}->{5}->{$_} || 0} @bases );
	    $size2ftype_count{$c}->{$type} = $sum;
	}
	close($rpt);	
    }
    open(my $rpt => ">$odir/$dbname.feature_sizes.dat") || die $!;
    print $rpt join("\t", qw(LENGTH), @kinds), "\n";
    for my $size ( sort { $a <=> $b } keys %size2ftype_count ) {
	print $rpt join("\t", $size, map { $size2ftype_count{$size}->{$_} || 0 } @kinds),"\n";       
    }
    close($rpt);
    print $R "pdf(\"$dbname\_reads_by_feature_and_size_combined.pdf\")\n";
    print $R "sizef <- read.table(\"$dbname.feature_sizes.dat\",header=T,sep=\"\\t\",row.names=1)\n";
    print $R "sizefp <- prop.table(as.matrix(sizef),margin=1)*100;\n";

    print $R "barplot(t(sizef),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname Reads mapped by size and feature type\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=sizef\$V1,legend=T,col=rainbow(",scalar @kinds,",start=0.1, end=.91),beside=F)\n";
    print $R "barplot(t(sizefp),xlab=\"Read Length\", ylab=\"Total # Reads\", main=\"$dbname Reads mapped by size (percent) and feature type\",space=0.1,cex.axis=0.8,las=1,cex=0.8,names=sizef\$V1,legend=T,col=rainbow(",scalar @kinds,",start=0.1, end=.91),beside=F)\n";    
}


sub compute_holes {
    my $data = shift;
    my $lens = shift;
    my $holes = {};
    while( my ($nm,$chrom) = each %{$data} ) {
	my $start = 1;
	for my $feat ( sort { $a->[0] <=> $b->[0] } @$chrom ) {
	    if( $feat->[0] > $start ) {
		push @{$holes->{$nm}}, [$start => $feat->[0] -1];
	    }
	    $start = $feat->[1]+1;
	}
	if( $start < $lens->{$nm} ) {
	    push @{$holes->{$nm}}, [$start => $lens->{$nm}];
	}
    }
    $holes;
}

