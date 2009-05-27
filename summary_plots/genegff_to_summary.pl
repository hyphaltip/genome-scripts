#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Bio::DB::SeqFeature;
my $debug = 0;

my $version = '2'; # 
my $method  = 'predicted';

GetOptions('v|verbose|debug!' => \$debug,
	   'version:i'        => \$version,
	   'm|method:s'       => \$method,
	   );

my $dir = shift;
my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'berkeleydb',
					 -dir     => $dir);

open(my $genefh => ">gene_summary.tab")  || die $!;
open(my $interfh => ">intergenic_summary.tab") || die $!;

print $genefh "#", join("\t", qw(gene transcript clone chrom strand pep_length 
			 spliced_length start stop chrom_start chrom_stop
			 genome_start genome_stop gmap source version method)),
    "\n";

print $interfh "#", join("\t", qw(chrom chrom_start chrom_stop length
				  left_gene left_strand 
				  right_gene right_strand)),"\n";

my %CHROMS;
{
    my $i = 0;
    for my $chrom ( sort { $a->id cmp $b->id } 
		    $db->get_features_by_type('scaffold') ) {
	$CHROMS{$chrom->seq_id} = [$i++,$chrom->length];
    }
}
my $offset = 0;
for my $chrom ( sort { $CHROMS{$a}->[0] <=> $CHROMS{$b}->[0] }
		keys %CHROMS ) {
    my $segment = $db->segment($chrom);
    my $lastgene;
    
    for my $gene ( sort { $a->start <=> $b->start } 
		   $segment->features('gene') ) {
	if( ! $lastgene ) {
	    if( $gene->start > 1 ) {
#		print $interfh join("\t", $gene->seq_id,
#				    1,$gene->start-1,
#				    $gene->start-1,
#				    'CHROM_START','+1', 
#				    $gene->id, 
#				    $gene->strand > 0 ? "+1" : "-1" ),"\n";
	    }
	} else {
	    my $distance = $gene->start-1 - $lastgene->end;
	    print $interfh join
		("\t", $gene->seq_id,
		 $lastgene->end+1,$gene->start-1,
		 $distance,
		 $lastgene->id,$lastgene->strand > 0 ? "+1" : "-1",
		 $gene->id,$gene->strand > 0 ? "+1" : "-1" ),"\n";
	}
	$lastgene = $gene;
	for my $mRNA ( $gene->get_SeqFeatures ) {
	    my $cds_len = 0;
	    for my $cds ( $mRNA->get_SeqFeatures('cds') ) {

		$cds_len += $cds->length;
	    }
	    print $genefh join("\t", $gene->id, $mRNA->id, 
			       $mRNA->seq_id, $mRNA->seq_id,
			       $gene->strand > 0 ? "+1" : "-1", ($cds_len/3) -1,
			       $cds_len,$gene->start, 
			       $gene->end, $gene->start,$gene->end, 
			       $offset + $gene->start,$offset + $gene->end, 
			       0, $gene->source_tag, $version, $method),"\n";
	}
	last if $debug;
    }
    last if $debug;
    $offset += $CHROMS{$chrom}->[1];
    if( $lastgene ) {
	if( $lastgene->end < $CHROMS{$chrom}->[1]) {
	    print $interfh join
		("\t", $lastgene->seq_id, 
		 $lastgene->end+1,$CHROMS{$chrom}->[1], 
		 ($CHROMS{$chrom}->[1] - $lastgene->end),
		 $lastgene->id, $lastgene->strand > 0 ? "+1" : "-1",
		 'CHROM_END','+1'),"\n";
	}
    }
}
