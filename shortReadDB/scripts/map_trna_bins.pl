#!/usr/bin/perl -w

=head1 NAME 

map_trna_bins - Summarize the trna data (need to run dump_feature_overlap 
    with the --tRNA option) (PREFIX_trna.dat)

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
use List::Util qw(max min sum);

my $bin_count = 10;
my $by_feature; # dump out a file for each tRNA feature instead of merged bins
my $debug = 0;
my $sep = "\t";
my $trimmed_only = 0;
my $aligned_length = 0;

GetOptions(
	   'v|verbose!'  => \$debug,
	   'b|bins:i'    => \$bin_count,
	   'f|feature!'  => \$by_feature,
	   'sep:s'            => \$sep,
	   'trimmed!'         => \$trimmed_only,
	   'aligned!'         => \$aligned_length,
	   );


my ($len_column,$seq_column,$rstrand_column,
    $fstrand_column,$matchstart_column,
    $matchend_column, $fstart_column, $fend_column,
    $fname_column) = qw(READ_TRIM_LENGTH
		        READ_SEQ 
			MATCH_STRAND
			FEATURE_STRAND
			MATCH_START
			MATCH_END
			FEATURE_START
			FEATURE_END
			FEATURE);

for my $file ( @ARGV ) {
    unless( $file =~ /(\S+)\.dat(\.gz)?$/ ) {
	warn("cannot grab stem from $file\n");
	next;
    }
    my ($stem) = $1;    
    my $in;
    if( $2 ) {
	open($in => "zcat $file |") || die $!; 
    } else {
	open($in => $file) || die $!;
    }
    my $i = 0;
    my %header = map { $_ => $i++ } split(/\s+/,<$in>);

    my $len_colidx        = $header{$len_column};
    my $seq_colidx        = $header{$seq_column};
    my $rstr_colidx       = $header{$rstrand_column};
    my $fstr_colidx       = $header{$fstrand_column};
    my $matchstart_colidx = $header{$matchstart_column};
    my $matchend_colidx   = $header{$matchend_column};
    my $featname_colidx   = $header{$fname_column};
    my $fstart_colidx     = $header{$fstart_column};
    my $fend_colidx       = $header{$fend_column};

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
	die("Cannot find column $fstrand_column\n");
    }
    unless( defined $featname_colidx ) {
	die("Cannot find column $fname_column\n");
    }
    unless( defined $fstart_colidx ) {
	die("Cannot find column $fstart_column\n");
    }
    unless( defined $fend_colidx ) {
	die("Cannot find column $fend_column\n");
    }
    
    my (%fs,@bins,%sizes);
    my %expression;
    $i = 0;
    while(<$in>) {
	my @row = split;
	my ($size,$seq) = ($row[$len_colidx],
			   $row[$seq_colidx]);
	if( $aligned_length ) {
	    $size = abs($row[$matchstart_colidx]-$row[$matchend_colidx]);
	}
	if( ! exists $fs{$row[$featname_colidx]} ) {
	    $fs{$row[$featname_colidx]}->{Opp}->{ abs($row[$fend_colidx] - $row[$fstart_colidx])} = [];
	    $fs{$row[$featname_colidx]}->{Same}->{ abs($row[$fend_colidx] - $row[$fstart_colidx])} = [];
	}
	$expression{$row[$featname_colidx]}->[$size]++;
#	$freqs{$size}{total}++;
#	$freqs{$size}->{'prime5'}->{substr($seq,0,1)}++;
#	$freqs{$size}->{'prime3'}->{substr($seq,$size-1,1)}++;
#	if( $row[$rstr_colidx] == $row[$fstr_colidx] ) {
#	    $freqs{$size}->{'SameStrand'}->{substr($seq,0,1)}++;
#	} else {
#	    $freqs{$size}->{'OppStrand'}->{substr($seq,0,1)}++;
#	}

	my $offset  = $row[$matchstart_colidx] - $row[$fstart_colidx];
	$offset = 0 if( $offset < 0 );
	#if( $by_feature ) {
	my $strandmatch = $row[$rstr_colidx] == $row[$fstr_colidx] ? 'Same' : 'Opp';
	for( my $j = 0;$j < $size;$j++ ) {
	    my $pos = $j+$row[$matchstart_colidx] - $row[$fstart_colidx];
	    $fs{$row[$featname_colidx]}->{$strandmatch}->{$pos}->[$size]++;
	}
	#} else {
	$bins[int($offset / $bin_count)]->{$strandmatch}->[$size]++;

	#}
	$sizes{$size}++;
	$i++;
	warn("$i records processed\n") if ( $debug && 
					    ( $i > 1_000 && ($i % 1_000) ==0));
	last if $i > 100_000 && $debug;
    }

    my @sizes = sort { $a <=> $b} keys %sizes;

    open(my $fh => ">$stem.tRNA_expressionn.tab") || die $!;
    print $fh join("\t", qw(tRNA TOTAL), map{ 'R_'. $_ } @sizes),"\n";
    for my $gene_set ( # schwartzian transformation
		   sort { $a->[0] <=> $b->[0] } # for speedy sorting on 
		   map { [sum( grep {defined} @{$expression{$_}}), $_ ] } # computed item
		   keys %expression ) {
	print $fh join("\t", $gene_set->[1], $gene_set->[0],
		       map { $expression{$gene_set->[1]}->[$_] || 0 } @sizes),"\n";
    }
    close($fh);
    warn("sizes are ",join(",", @sizes),"\n") if $debug;
    $stem = sprintf("%s_%s",$aligned_length ? 'aln' : 'trim',$stem);
    open(my $R_fh => ">$stem.tRNA_bins.R") || die $!;
    # tRNA overall counts, normalized to bins
    printf $R_fh 'pdf("%s.tRNA_bins_counts.pdf")'."\n",$stem;    
    printf $R_fh 'counts <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
    "$stem.tRNA_bins_counts_SS.tab";
# same and opp strand feature breakdown
    printf $R_fh 'barplot(t(counts),xlab="BIN", ylab="# of reads", main="%s Read Dist on all tRNA features (Same Strand)",space=0.1,legend.text=T,beside=F,col=rainbow(%d,start=0.1,end=0.91),cex.main=0.8,cex.names=0.7)'."\n",$stem,scalar @sizes;
    printf $R_fh 'counts <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
    "$stem.tRNA_bins_counts_OS.tab";
    
    printf $R_fh 'barplot(t(counts),xlab="BIN", ylab="# of reads", main="%s Read Dist on all tRNA features (Opposite Strand)",space=0.1,legend.text=T,beside=F,col=rainbow(%d,start=0.1,end=0.91),cex.main=0.8,cex.names=0.7)'."\n",$stem, scalar @sizes;
        
    open(my $fh_ss => ">$stem.tRNA_bins_counts_SS.tab") || die $!;
    print $fh_ss join("\t", qw(BIN),@sizes ), "\n";
    open(my $fh_os => ">$stem.tRNA_bins_counts_OS.tab") || die $!;
    print $fh_os join("\t", qw(BIN),@sizes ), "\n";
    my $binid = 0;
    for my $bin ( @bins ) {
	print $fh_ss join("\t", $binid, 
			  map { $bin->{Same}->[$_] || 0} @sizes),"\n";
	print $fh_os join("\t", $binid, 
			  map { $bin->{Opp}->[$_] || 0} @sizes),"\n";
	$binid++;
    }
    close($fh_ss);
    close($fh_os);

    printf $R_fh 'pdf("%s.tRNA_density.pdf")'."\n",$stem;
    
    for my $fname ( sort keys %fs ) {
	for my $strandmatch ( qw(Same Opp) ) {
	    open(my $fh => ">$stem.tRNA_density_$strandmatch.$fname.tab") || die $!;
	    print $fh join("\t", qw(BASE ),@sizes),"\n";

	    printf $R_fh 'density <-read.table("%s",sep="\t",header=T,row.names=1)'."\n",
	    "$stem.tRNA_density_$strandmatch.$fname.tab";

	    printf $R_fh 'barplot(t(density),xlab="tRNA position (nt)", ylab="# of reads at position", main="%s read Density, %s Strand",space=0.1,legend.text=T,beside=F,col=rainbow(%d,start=0.1,end=0.91),cex.main=0.8,cex.names=0.7);'."\n",
			$fname,$strandmatch,scalar  @sizes;

	    for my $pos ( sort { $a <=> $b } keys %{$fs{$fname}->{$strandmatch}} ){
		print $fh join("\t", $pos,
			       map { $fs{$fname}->{$strandmatch}->{$pos}->[$_] || 0 }
			       @sizes),"\n";
	    }    
	    close($fh);
	}
    }
    close($R_fh);
}

