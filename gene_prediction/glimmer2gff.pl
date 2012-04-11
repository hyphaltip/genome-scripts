#!/usr/bin/perl -w
use strict;
use Bio::Tools::Glimmer;
use Bio::Tools::GFF;
use Bio::DB::Fasta;
use IO::String;

my $genome = shift;
my $db = Bio::DB::Fasta->new($genome);
my $dir = shift;
opendir(DIR,$dir) || die $!;
for my $file ( readdir(DIR) ) {
    if( $file =~ /(\S+)\.predict/) {
	my $id = $1;
	my $len = $db->length($id);
	if( ! $len ) {
	    warn("no length for $id\n");
	}
	my $out = Bio::Tools::GFF->new(-gff_version => 3,
				       -file => ">$id.gff3");

# file
	my $parser = Bio::Tools::Glimmer->new(-format    => 'Glimmer',
					      -seqlength => $len,
					      -seqname   => $id,
					      -file      => "$dir/$file" );
	while(my $gene = $parser->next_prediction()) {
	    $out->write_feature($gene);	
	}
    }
}
