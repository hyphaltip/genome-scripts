#!/usr/bin/perl -w
use strict;
my $out = Bio::Tools::GFF->new(-gff_version => 3);
for my $file ( @ARGV) {

my $gw = Bio::Tools::Genewise->new(-fh => $file);
    while (my $gene = $gw->next_prediction){
	$gene->primary_tag('gene');
	$gene->source_tag('Genewise');	
	$out->write_feature($gene);
	my @transcripts = $gene->transcripts;
	foreach my $t(@transcripts){
	    $t->primary_tag('mRNA');
	    $t->source_tag('Genewise');
	    $out->write_feature($t);
	    my @exons =  $t->exons;	    
	    foreach my $e(@exons){
		$e->remove_tag('supporting_feature');
		$e->source_tag('Genewise');
		$e->primary_tag('CDS');
		$out->write_feature($e);
	    }
	}
    }
}
