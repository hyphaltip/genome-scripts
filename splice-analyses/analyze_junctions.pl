#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::SeqFeature::Store;
my ($genome_junctions,@junctions) = @ARGV;

my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'berkeleydb',
					  -temp    => 1);

open(my $fh => $genome_junctions) || die $!;
my $n = 0;
while(<$fh>) {
    chomp;
    my ($chrom,$start,$end,$strand) = split(/\t/,$_);
    
    $dbh->add_SeqFeature( $dbh->new_feature(-seq_id => $chrom,
					    -start  => $start,
					    -end    => $end,
					    -strand => $strand eq '+' ? 1 : -1,
					    -primary_tag => 'junction',
					    -source_tag  => 'observed'));
    $n++;
}

for my $tophat_junctions ( @junctions ) {
    open($fh => $tophat_junctions) || die $!;
    open(my $ofh => ">$tophat_junctions.nonoverlap.bed") || die $!;
    my $i =0;
    while(<$fh>) {
	next if /^track/ || /^\s+$/;
	chomp;
	my ($chrom,$start,$end,$name,$depth,$strand) = split(/\t/,$_);
	if( ++$i % 1000 == 0 ) {
	    warn("$i junctions processed\n");
	}
	next if( $dbh->features(-type=> 'junction:observed',
				-start => $start,
				-end   => $end,
				-seq_id=>$chrom,
				-range_type => 'overlaps') );
	print $ofh join("\t",$chrom,$start,$end,$name,$depth,$strand), "\n";
    }
}

