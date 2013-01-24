#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::Perl qw(revcom_as_string);

my $window = 50;

my $trimseq = shift || 'Ccin_hyphaltip_RNAseq.trinity.trimpoly.seqs';
my $alnfolder = shift || 'alns';
my $genome = shift || 'coprinopsis_cinerea_okayama7_genome.fasta';
my $db = Bio::DB::Fasta->new($genome);

my $in = Bio::SeqIO->new(-format => 'fasta',
			 -file   => $trimseq);
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -file   => ">exonerate_based_polyA.fas");

my $out_neg = Bio::SeqIO->new(-format => 'fasta',
			  -file   => ">exonerate_based_polyA_negative.fas");
while( my $seq = $in->next_seq ) {
    my $id = $seq->display_id;
    my $desc = $seq->desc;
    if( -f "$alnfolder/$id.exonerate" ) {
	my $searchio = Bio::SearchIO->new(-format => 'exonerate',
					  -file   => "$alnfolder/$id.exonerate");
	while( my $r = $searchio->next_result ) {
	    while( my $h = $r->next_hit ) {
		print join("\t",$r->query_name, 
			   $desc,
			   $h->start('query'),
			   $h->end('query'),
			   $h->name,
			   $h->start('hit'),
			   $h->end('hit')),"\n";

		if( $desc eq '5prime' ) {
		    my $target_window = $db->seq($h->name,
						 $h->start('hit') - $window => 
						 $h->start('hit') + $window -1);
		    next if ( $target_window =~ /T{5,}/ ||
			      $target_window =~ /A{5,}/ );
		    
		    $out->write_seq(Bio::Seq->new(-seq => $target_window,
						  -id  => $id,
						  -desc => sprintf("%s:%s",
								   $h->name,
								   $h->start('hit'))));
		    
		    my $up = $db->seq($h->name,
				      $h->start('hit') - 4 * $window => 
				      ($h->start('hit') - 3 * $window - 1) );
		    
		    $out_neg->write_seq(Bio::Seq->new(-id => "UP_".$id,
						      -seq => $up));

		    $up = $db->seq($h->name,
				   $h->start('hit') + 3 * $window => 
				   ($h->start('hit') + 4 * $window -1) );
		    
		    $out_neg->write_seq(Bio::Seq->new(-id => "UP2_".$id,
						      -seq => $up));

		} else {
		    my $target_window = $db->seq($h->name,
						 $h->end('hit') - $window =>
						 $h->end('hit') + $window -1);
		    next if ( $target_window =~ /T{5,}/ ||
			      $target_window =~ /A{5,}/ );
		    
		    $out->write_seq(Bio::Seq->new(-seq => $target_window,
						  -id  => $id,
						  -desc => sprintf("%s:%s",
								   $h->name,
								   $h->end('hit'))));
		    my $down = $db->seq($h->name,
					$h->end('hit') + 3*$window => 
					$h->end('hit') + 4*$window -1);
		    $out_neg->write_seq(Bio::Seq->new(-id => "DN_".$id,
						      -seq => $down));
		    $down = $db->seq($h->name,
				     $h->end('hit') - 4*$window => 
				     $h->end('hit') - 3*$window -1);
		    $out_neg->write_seq(Bio::Seq->new(-id => "DN2_".$id,
						      -seq => $down));
		}
	    }
	}    
    }
}

