#!/usr/bin/perl -w
use strict;
use Bio::Tools::tRNAscanSE;
use Bio::Tools::GFF;

my $out = Bio::Tools::GFF->new(-gff_version => 3);
my $parser = Bio::Tools::tRNAscanSE->new(-file => shift);

# parse the results
while( my $gene = $parser->next_prediction ) {
    my ($id) = $gene->get_tag_values('ID');
    $gene->primary_tag('gene');
    $out->write_feature($gene);
    my $i = 0;
     for my $trans ( $gene->get_SeqFeatures() ) {	
	 $trans->primary_tag('exon');
	 $out->write_feature($trans,$trans->get_SeqFeatures());
	 $i++;
    }

}

