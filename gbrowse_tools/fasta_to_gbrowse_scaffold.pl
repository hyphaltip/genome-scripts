#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
use Getopt::Long;
my $version = 3;

GetOptions('v|version:i' => \$version);

my $in = new Bio::SeqIO(-fh => \*ARGV,
			-format=> 'fasta');
my $out = new Bio::Tools::GFF(-gff_version=>$version);
print "##gff-version $version\n" if $version == 3;
while( my $seq = $in->next_seq ) {
	my $group;
     if( $version == 3 ){
	$group = sprintf("ID=%s;Name=%s", $seq->display_id, $seq->display_id);
     } else {
	$group = sprintf("Scaffold %s", $seq->display_id);
     }

   print join("\t", $seq->display_id,
   		    'chromosome',
		    'scaffold',
		    1,
		    $seq->length,
		    '.',
		    '.',
		    '.',
		    $group,
		    ),"\n";
}
