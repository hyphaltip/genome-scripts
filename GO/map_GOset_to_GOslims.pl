#!/usr/bin/perl -w
use strict;

=head1 NAME

map_GOset_to_GOslim.pl - map a datafile of GO terms into the appropriate GO slim categorizes

=head1 USAGE

map_GOset_to_GOslim.pl -i organism.go -db gene_ontology.obo > organism.go_slim

=head1 SYNOPSIS

This script will take a datafile of GO terms and map them into GO slim terms for a defined GO slim to effectively reduce the complexity of the terms in the file.

It currently expects the simple GeneSetEnrichment table format of

  GO:12345\tIEA\tGENE_NAME

Where GO:12345 is any GO term, IEA is the attribution of annotation, could be any of the curated or automated descriptions, and the GENE_NAME is simple string with the gene identifier or accession number.

--go or --db or --godb provides gene_ontology.obo file
-i -or --input         provides the GO to gene table format described above

=head1 AUTHOR

Jason Stajich E<lt>jason.stajich[at]ucr.eduE<gt>

=cut

use Getopt::Long;
# GO module
use GO::Parser;

my $godbfile = '/srv/projects/db/GO/current/gene_ontology.obo';
my $GOslim_set = 'goslim_yeast';
my $input;
GetOptions
    ('go|db|godb:s' => \$godbfile,
     'i|input:s'    => \$input,
     'slim:s'       => \$GOslim_set,
     'h|help'      => sub { 
	 exec('perldoc',$0);
	 exit;
     }
    );
my $fh;

if( ! defined $input || ! -f $input ) {
    open($fh => \*ARGV) || die "cannot open input from ARGV\n";
} else {
    open($fh => $input) || die "cannot open $input: $!";
}

my $parser = new GO::Parser({handler=>'obj', use_cache => 1}); # create parser object
$parser->parse($godbfile); # parse the database

my $graph = $parser->handler->graph;  # get GO::Model::Graph object

while(<$fh>) {
    chomp;
    my ($interm,$src,$gene) = split(/\t/,$_);
    my $term = $graph->get_term($interm);
    my %parent_terms = &gather_parents_slim($term,$graph);
    for my $acc ( sort keys %parent_terms ) {
	my $i_term = $parent_terms{$acc};
	print join("\t", $acc, $src, $gene, $i_term->name), "\n";    
    }
}

# recursively gather parents, checking if they are in the right
# go slim, stop as soon as we see a parent with this

sub gather_parents_slim {
    my ($term,$p_graph) = @_;
    my %terms;
    if( $term->in_subset($GOslim_set) ) {
	%terms = ($term->acc => $term );
    } else {
	my $rels = $p_graph->get_parent_relationships($term->acc);
	my @terms_here;
	for my $rel ( @$rels ) {
	    next unless $rel->type eq 'is_a';
	    #warn( join(" ", $rel->subject_acc, $rel->type, $rel->object_acc),"\n");	    
	    my $object = $p_graph->get_term($rel->object_acc);
	    push @terms_here, $object;
	    if(  $object->in_subset($GOslim_set) && 
		 ! exists $terms{$object->acc} ) {		
		$terms{$object->acc} = $object;
	    }
	}
	# so if any parents of this term (eg made it into @terms) 
	# then we are done looking for a parent term that is in a GO-slim
	if( ! keys %terms ) {
	    # if not we want to get the parent terms of all of these
	    # and find one that is in GO slim
	    for my $t ( @terms_here ) {
		my %r = &gather_parents_slim($t,$p_graph);
		# combine the returned values from different
		# result sets
		
		while( my ($r_k,$r_v) = each %r ) {
		    $terms{$r_k} = $r_v;
		}
	    }
	}
    }
    %terms;
}
