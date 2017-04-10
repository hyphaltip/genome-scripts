#!/usr/bin/perl 
use strict;

use Bio::TreeIO;

my $infile = shift;
my $out = Bio::TreeIO->new(-format => 'newick', 
			   -internal_node_id => 'bootstrap');

my $in = Bio::TreeIO->new(-format => 'newick', -file => $infile);

while(my $t = $in->next_tree ) {
    $out->write_tree($t);
#    for my $node ( $t->get_nodes ) {
#	printf "id: %s bootstrap: %s\n", 
#	$node->id || '', 
#	$node->bootstrap || '', 
#	"\n";
#    }
}
