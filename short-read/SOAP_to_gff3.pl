#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(sum);
my $debug = 0;
my $manyhits = 10;
my $type   = 'SOAP';
my $source = 'Solexa_miRNA';
my $condtype = 'SOAP_condensed';

GetOptions(
	   't|type:s'         => \$type,
	   's|source:s'       => \$source,
	   'c|condense:s'     => \$condtype,
	   'v|verbose|debug!' => \$debug,
	   );
my %compress = ('bz2' => 'bzcat',
		'gz'  => 'zcat');

my @ext = qw(multihit fewhit condensed);

for my $file ( @ARGV ) {
    my ($fh,$ofh);
    my $stem = $file;
    if( $file =~ /(.+)\.(bz2|gz)$/ ) {
	open($fh, "$compress{$2} $file |") || die $!;	
	$stem = $1;
    } else {
	open($fh, "< $file") || die $!;
    }
    $stem =~ s/\.out$//;
    my %fh;
    for my $ex ( @ext ) {
	open($fh{$ex},">$stem.$ex.gff3") || die "$stem.$ex.gff3: $!";
	$ofh = $fh{$ex};
	printf $ofh "##gff-version 3\n##date %s\n","".localtime(time);
    }

    for my $dir ( qw(first term) ) {
	open($fh{"$dir.sum.low"},">$stem.$dir.simple_summary.tab") || die "$stem.simple_summary: $!";
	open($fh{"$dir.sum.all"},">$stem.$dir.summary.tab") || die "$stem.simple_summary: $!";
    }
    my $i =0;
    my %locs;
    my (%bases, %lowcopybases);
    while(<$fh>) {
	chomp;
	my ($read,$seq,$qual, $no_hits,
	    undef, $trimmed_len,
	    $strand, $chrom, $start,
	    $mm,@mms) = split(/\t/,$_);
	my $firstbase = uc substr($seq,0,1);
	my $lastbase = uc substr($seq,-1,1);
	
	$bases{'first'}->{$firstbase}->{$trimmed_len}++;
	$bases{'term'}->{$lastbase}->{$trimmed_len}++;
	if( $no_hits <= $manyhits ) {
	    $lowcopybases{'first'}->{$firstbase}->{$trimmed_len}++;
	    $lowcopybases{'term'}->{$lastbase}->{$trimmed_len}++;
	}
	my $end = $start + $trimmed_len;
	my $locstr = sprintf("%d:%d:$strand",sort { $a <=> $b } $start,$end);
	$locs{$chrom}->{$locstr}->{'reads'}++;
	$locs{$chrom}->{$locstr}->{'hits'} += $no_hits;
	my $mmstr = '';
	if( $mm ) {
	    for ( @mms ) { s/\->/:/g }	
	    $mmstr = sprintf(";NumMismatches=%d;Mismatches=\"%s\"",
			     $mm, join(",",@mms));
	}
	my $ttype;
	if( $no_hits > $manyhits ) {
	    $ofh = $fh{'multihit'};
	    $ttype = $type. "_multi";
	} else {
	    $ofh = $fh{'fewhit'};
	    $ttype = $type. "_few";
	}
	print $ofh join("\t",
		   $chrom, $source,$ttype,			
		   $start, $end,
		   $trimmed_len,
		   $strand, '.',
		   sprintf("ID=%s;LeadingBase=%s;NumHits=%d%s",
			   $read,$firstbase,$no_hits,$mmstr)),"\n";
	last if $debug && $i++ > 10;
    }
    close($fh);
    $ofh = $fh{'condensed'};
    my $c = 0;
    for my $chr ( sort keys %locs ) {
	for my $pos ( sort { $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2] } 
		      map { [$_,split(/:/,$_)] } 
		      keys %{$locs{$chr}} ) {
	    my $d = $locs{$chr}->{$pos->[0]};
	    next unless $d->{'reads'} > 1; # skip singletons?
	    print $ofh join("\t",$chr, $source,$condtype,
			    $pos->[1],$pos->[2],
			    $d->{'reads'}, # score is # of reads
			    $pos->[3], # strand
			    '.', # frame is empty
			    sprintf("ID=sRNAClus%05d;GenomicHits=%d",$c++,
				    $d->{'hits'})),"\n";
	
	}
    }

    for my $m ( qw(low all) ) {
	for my $dir ( qw( first term) ) {
	    $ofh = $fh{"$dir.sum.$m"};
	    my @bases;
	    my %freqs;
	    for my $base ( sort keys %{$bases{$dir}} ) {
		next if $base !~ /[CAGT]/i;
		while( my ($len,$count) = each %{$bases{$dir}->{$base}} ) { 	    
		    $freqs{$len}->{$base} = $count;
		}
		push @bases, $base;
	    }
	    
	    print $ofh join("\t", 'LEN', @bases),"\n";
	    for my $len ( sort { $a <=> $b } keys %freqs ) {
		my $sum = sum ( values %{$freqs{$len}} );
		print $ofh join("\t", $len, map 
			   { sprintf("%.6f",
				     $freqs{$len}->{$_}/$sum) }
			   @bases),"\n";
	    }
	}
    }
}
