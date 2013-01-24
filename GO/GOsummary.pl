#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);

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

my %gather;
my %lookup;
my %total_genes;
while(<$fh>) {
    chomp;
    my ($interm,$src,$gene) = split(/\t/,$_);
    my $term = $graph->get_term($interm);
    $lookup{$term->acc} = $term->name unless exists $lookup{$term->acc}; 
    $gather{$term->acc}->{$gene}++;
    $total_genes{$gene}++;
}

my $total_gene_count = scalar %total_genes;
print join("\t", qw(TERM DESC COUNT PERCENTAGE)),"\n";
for my $term ( sort { $gather{$b} <=> $gather{$a} }
	       keys %gather ) {
    my $gene_count = scalar keys %{$gather{$term}};
    print join("\t", $term, $lookup{$term}, $gene_count, 
	       sprintf("%.2f",100 * ($gene_count / $total_gene_count))), "\n";
}
