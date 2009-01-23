#!/usr/bin/perl -w

=head1 NAME

feat2sum - Summarize the size class distribution of short Reads
sub-divided by the types of features they overlap.

=head1 DESCRIPTION

To run this script
 cd  taylor_projects/sReadDB/results
 perl ../scripts/feat2sum.pl *.dat
 for a in `ls *.R`
 do
  R --no-save < $a
 done


Will generate all the PDFs.

=cut

use strict;
use Getopt::Long;
use Env qw(HOME);

my ($debug);

my ($len_column,$seq_column,$rstrand_column,
    $fstrand_column,$matchstart_column,
    $matchend_column) = qw(READ_TRIM_LENGTH
			   READ_SEQ 
			   MATCH_STRAND
			   FEATURE_STRAND
			   MATCH_START
			   MATCH_END);
my $sep = "\t";
my $trimmed_only = 0;
my $aligned_length = 0;
GetOptions(
	   'v|verbose|debug!' => \$debug,
	   'sep:s'            => \$sep,
	   'trimmed!'         => \$trimmed_only,
	   'aligned!'         => \$aligned_length,
	   );

my @nt = qw(A C G T);
for my $file ( @ARGV ) {
    
    unless( $file =~ /(\S+)\.dat(\.gz)?/ ) {
	warn("cannot grab stem from $file\n");
	next;
    }
    my ($stem,$ext) = ($1,$2);    
    my $in;
    if( $ext ) { # handle compressed files
      open($in => "zcat $file |") || die $!;
    } else {
      open($in => $file) || die $!;
    }
    my $i = 0; 
    my %header = map { $_ => $i++ } split(/\s+/,<$in>);
    my $len_colidx = $header{$len_column};
    my $seq_colidx = $header{$seq_column};
    my $rstr_colidx = $header{$rstrand_column};
    my $fstr_colidx = $header{$fstrand_column};
    my $matchstart_colidx = $header{$matchstart_column};
    my $matchend_colidx   = $header{$matchend_column};

    unless( defined $len_colidx ) { 
	die("Cannot find column $len_column\n");
    }
    unless( defined $seq_colidx ) { 
	die("Cannot find column $seq_column\n");
    }
    unless( defined $rstr_colidx ) {
	die("Cannot find column $rstrand_column\n");
    }
    unless( defined $fstr_colidx ) {
	warn("Cannot find column $fstrand_column\n");
    }
    
    my %freqs;
    $i = 0;
    while(<$in>) {
	my @row = split;
	my ($size,$seq) = ($row[$len_colidx],
			   $row[$seq_colidx]);
	if( $aligned_length ) {
	    $size = abs($row[$matchstart_colidx]-$row[$matchend_colidx]);
	}
	next unless defined $size;
	$freqs{$size}{total}++;
	$freqs{$size}->{'prime5'}->{substr($seq,0,1)}++;
	$freqs{$size}->{'prime3'}->{substr($seq,$size-1,1)}++;
	if( defined $fstr_colidx )  {
	 if( $row[$rstr_colidx] == $row[$fstr_colidx] ) {
	     $freqs{$size}->{'SameStrand'}->{substr($seq,0,1)}++;
	 } else {
	     $freqs{$size}->{'OppStrand'}->{substr($seq,0,1)}++;
	 }
	}
	$i++;
	warn("$i records processed\n") if ( $debug && 
					    ( $i > 1_000 && ($i % 1_000) ==0));
	last if $i > 10_000 && $debug;
    }
    $stem = sprintf("%s_%s",$aligned_length ? 'aln' : 'trim',$stem);
    open(my $R_fh => ">$stem.R") || die $!;

    printf $R_fh 'pdf("%s.counts.pdf")'."\n",$stem;    
    printf $R_fh 'counts <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
    "$stem.counts.tab";
    printf $R_fh 'barplot(t(counts),xlab="Read Length", ylab="# of reads", main="%s Size Distribution",space=0.1,legend.text=T,beside=F,col=rainbow(5,start=0.1,end=0.91),cex.main=0.8,cex.names=0.9)'."\n",$stem;
    
    open(my $fh => ">$stem.counts.tab") || die $!;
    print $fh join("\t", qw(SIZE),@nt), "\n";
    for my $size ( sort { $a <=> $b } keys %freqs ) {
#	print $fh join("\t", $size, $freqs{$size}->{total}),"\n";
	print $fh join("\t", $size, map { $freqs{$size}->{prime5}->{$_} || 0} @nt),"\n";
    }
    printf $R_fh 'pdf("%s.5prime_freq.pdf")'."\n",$stem;    
    printf $R_fh 'prime5 <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
    "$stem.5prime_freq.tab";
    printf $R_fh 'barplot(t(prime5),xlab="Read Length", ylab="%% of reads", main="%s 5\' Base Freq",space=0.1,legend.text=T,beside=F,col=rainbow(5,start=0.1,end=0.91),cex.main=0.8,cex.names=0.9)'."\n",$stem;

    open($fh => ">$stem.5prime_freq.tab") || die $!;
    print $fh join("\t", qw(SIZE), @nt), "\n";
    for my $size ( sort { $a <=> $b } keys %freqs ) {
	print $fh join("\t", $size, map { 
	    sprintf("%.3f",100 * ($freqs{$size}->{prime5}->{$_}||0) /
		    $freqs{$size}->{total})
	    } @nt),"\n";
    }
    close($fh);

    printf $R_fh 'pdf("%s.3prime_freq.pdf")'."\n",$stem;    
    printf $R_fh 'prime3 <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
    "$stem.3prime_freq.tab";
    printf $R_fh 'barplot(t(prime3),xlab="Read Length", ylab="%% of reads", main="%s 3\' Base Freq",space=0.1,legend.text=T,beside=F,col=rainbow(5,start=0.1,end=0.91),cex.main=0.8,cex.names=0.9)'."\n",$stem;

    open($fh => ">$stem.3prime_freq.tab") || die $!;
    print $fh join("\t", qw(SIZE), @nt), "\n";
    for my $size ( sort { $a <=> $b } keys %freqs ) {
	print $fh join("\t", $size, map { 
	    sprintf("%.3f",100 * ($freqs{$size}->{prime3}->{$_} || 0)/
		    $freqs{$size}->{total}) 
	    } @nt),"\n";
    }
    close($fh);
    if( defined $fstr_colidx) {
     printf $R_fh 'pdf("%s.StrandBias.pdf")'."\n",$stem;    
     printf $R_fh 'sstrand <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
     "$stem.samestrand_counts.tab";
     printf $R_fh 'barplot(t(sstrand),xlab="Read Length", ylab="# of reads", main="%s 5\' Base Freq, same strand as feature",space=0.1,legend.text=T,beside=F,col=rainbow(5,start=0.1,end=0.91),cex.main=0.8,cex.names=0.9)'."\n",
    $stem;
    printf $R_fh 'ostrand <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
    "$stem.oppstrand_counts.tab";
    printf $R_fh 'barplot(t(ostrand),xlab="Read Length", ylab="# of reads", main="%s 5\' Base Freq, OppStrand as feature",space=0.1,legend.text=T,beside=F,col=rainbow(5,start=0.1,end=0.91),cex.main=0.8,cex.names=0.9)'."\n",
    $stem;
    
    open($fh => ">$stem.samestrand_counts.tab") || die $!;
    print $fh join("\t", qw(SIZE),@nt), "\n";
    for my $size ( sort { $a <=> $b } keys %freqs ) {
#	print $fh join("\t", $size, $freqs{$size}->{total}),"\n";
	print $fh join("\t", $size, map { $freqs{$size}->{SameStrand}->{$_} || 0} @nt),"\n";
    }
    close($fh);

    open($fh => ">$stem.oppstrand_counts.tab") || die $!;
    print $fh join("\t", qw(SIZE),@nt), "\n";
    for my $size ( sort { $a <=> $b } keys %freqs ) {
#	print $fh join("\t", $size, $freqs{$size}->{total}),"\n";
	print $fh join("\t", $size, map { $freqs{$size}->{OppStrand}->{$_} || 0} @nt),"\n";
    }
    close($fh);
   }	
}
