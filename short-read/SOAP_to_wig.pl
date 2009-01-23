#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
my $debug = 0;
my $genome;
my $name = 'Solexa';
my $desc = 'SOAP mapped Solexa reads';
GetOptions(
	   'v|verbose|debug!' => \$debug,
	   'db|g|genome:s'       => \$genome,
	   'n|name:s'         => \$name,
           'd|desc:s'         => \$desc,
	   );
my %compress = ('bz2' => 'bzcat',
		'gz'  => 'zcat');

for my $file ( @ARGV ) {
    my ($fh,$ofh);
    my $stem = $file;
    if( $file =~ /(.+)\.(bz2|gz)$/ ) {
	open($fh, "$compress{$2} $file |") || die $!;	
	$stem = $1;
    } else {
	open($fh, "< $file") || die $!;
    }

    my $g = Bio::SeqIO->new(-format => 'fasta',
			    -file   => $genome);
    my %chroms;
    while( my $s = $g->next_seq ) {
	$chroms{$s->display_id}->[$s->length] = 0;
    }
    $stem =~ s/\.out$//;
    while(<$fh>) {
	chomp;
	my ($read,$seq,$qual, $no_hits,
	    undef, $trimmed_len,
	    $strand, $chrom, $start,
	    $mm,@mms) = split(/\t/,$_);
	for my $b ( $start..($start+$trimmed_len -1) ) {
	    $chroms{$chrom}->[$b]++;
	}
    }    
    close($fh);    
    open($ofh, ">$stem.wig") || die $!;
    printf $ofh "track type=wiggle_0 name=\"%s\" description=\"%s\"\n",$name,$desc;
    for my $chrom ( sort keys %chroms ) {
	print $ofh "fixedStep chrom=$chrom start=1 step=1\n";
	my $ar = $chroms{$chrom};
	my $i = 0;
	for ( @{$ar} ) {
	    if( $i ) {
		printf $ofh "%d\n",($_ || 0);
	    }
	    $i=1; # skip first since it is 0
	}
    }
    close($ofh);
}
