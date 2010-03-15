#!/usr/bin/perl -w
use strict;
use File::Path 'rmtree','mkpath';
use Bio::DB::SeqFeature::Store;

my $dS = 'plot/coprinus_avg_ds.tab';
my $recombo = 'plot/recombination_rates.tab';
rmtree('/scratch/dS');
mkpath('/scratch/dS');

my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor=> 'berkeleydb',
					  -dsn    => '/scratch/dS/',
					  -create => 1);
open(my $fh => $recombo ) || die $!;
my @hdr = split(/\s+/,<$fh>);
my %recombo_rates;

while(<$fh>) {
    my ($chrom,$start,$stop,$type,$tables,$designation) = split;
#    warn("chrom is $chrom, start is $start, end is $stop\n");
    $dbh->new_feature(-seq_id      => $chrom,
		      -start       => $start,
		      -end         => $stop,
		      -primary_tag => $designation,
		      -source      => 'recombination');

#    push @{$recombo_rates{$chrom}}, [$start,$stop,$designation];
}
close($fh);

open($fh => $dS) || die $!;
@hdr = split(/\s+/,<$fh>);
my %dS;

while(<$fh>) {
    my ($model,$count,$dS,$dN,$chrom,$start) = split;
    my @f = $dbh->get_features_by_location(-seq_id => $chrom,
					   -start  => $start,
					   -end    => $start+100);
    for my $f ( @f ) {
	my $type = $f->primary_tag;
	push @{$dS{$type}}, $dS;
    }
}
close($fh);
while( my ($type,$vals) = each %dS) {
    open(my $ofh => ">dS_$type.tab") || die $!;
    print $ofh join("\n", @$vals), "\n";
}
