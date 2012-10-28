#!/usr/bin/perl -w
use strict;
use File::Spec;
use Env;

use strict;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OntologyParser;
use GO::OntologyProvider::OboParser;
use GO::TermFinderReport::Text;

use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);

use Getopt::Long;
$|=1;

my $godir = "$HOME/lib/GO";
my $totalNum = 10_000;
my $assocfile;
my $cutoff = 0.05;

my $usage = "list_test_GOenrich: -a associationfile [-go godir] [-n totalnumgenes]\n";

GetOptions(
	   'a|assoc:s'    => \$assocfile,
	   'go|godir:s'   => \$godir,
	   'n|totalnum:s' => \$totalNum,
	   'c|cutoff:s'   => \$cutoff,
	   );

die("no Association file provided!\n$usage\n") unless( defined $assocfile && 
						       -f $assocfile);

my $process   = GO::OntologyProvider::OboParser->new(ontologyFile => "$godir/gene_ontology.1_2.obo",  aspect       => 'P');
my $component = GO::OntologyProvider::OboParser->new(ontologyFile => "$godir/gene_ontology.1_2.obo",  aspect	=> 'C');
my $function  = GO::OntologyProvider::OboParser->new(ontologyFile => "$godir/gene_ontology.1_2.obo", aspect  => 'F');

my @files = @ARGV;

# now get our annotation file and number of genes
# now set up the objects we need


my $annotation = GO::AnnotationProvider::AnnotationParser->new('annotationFile' => $assocfile);

my $report = GO::TermFinderReport::Text->new();

for my $file ( @files ) {
    warn( "Analyzing $file\n");

    my @genes;
    open(IN, $file) || die $!;
    while(<IN>) {
	my ($gene) = split;
	push @genes, $gene;
    }
    close(IN);
    my %termFinder;
    $termFinder{'P'} = GO::TermFinder->new(annotationProvider=> $annotation,
					   ontologyProvider  => $process,
					   totalNumGenes     => $totalNum,
					   aspect            => 'P');
	
    
    $termFinder{'C'} = GO::TermFinder->new(annotationProvider=> $annotation,
					   ontologyProvider  => $component,
					   totalNumGenes     => $totalNum,
					   aspect            => 'C');
	
    $termFinder{'F'} = GO::TermFinder->new(annotationProvider=> $annotation,
					   ontologyProvider  => $function,
					   totalNumGenes     => $totalNum,
					   aspect            => 'F');

    open(OUT, ">$file.GO_report") || die "Cannot open $file.GO_report: $!";
    open(WARN, ">$file.warn") || die "Cannot open $file.warn: $!";
    my (@list, @notFound, @ambiguous);    
    CategorizeGenes(annotation  => $annotation,
		    genes       => \@genes,
		    ambiguous   => \@ambiguous,
		    unambiguous => \@list,
		    notFound    => \@notFound);
    
    if (@list) {
	print WARN (scalar @list), " gene(s) will be considered:\n";

	foreach my $gene (@list){
	    print OUT $gene, "\t", $annotation->standardNameByName($gene), "\n";		
	}
	print OUT "\n";
    } else{
	print WARN "None of the gene names were recognized\n";
	print WARN "They were:\n\n";
	print WARN join("\n", @notFound), "\n";
	# close WARN;
	next;
    }
    
    if (@ambiguous) {
	print WARN scalar @ambiguous, " gene(s) are ambiguously named and unused\n";
	
	print WARN join("\n", @ambiguous), "\n\n";
    }
    if (@notFound) {
	print WARN "The following gene(s) were not recognized, and will not be considered:\n\n";
	print WARN join("\n", @notFound), "\n\n";
    }
    my %go_seen;
    for my $termFinder (map { $termFinder{$_} } qw(P C F) ) {
	print OUT "Finding terms for ", $termFinder->aspect, "\n\n";
	my @pvalues = $termFinder->findTerms(genes        => \@list,
					     calculateFDR => 1);
	
	my $numHypotheses = $report->print(pvalues  => \@pvalues,
					   numGenes => scalar(@list),
					   totalNum => $totalNum,
					   cutoff   => $cutoff,
					   fh       => \*OUT);
	

	# if they had no significant P-values
	if ($numHypotheses == 0){
	    print WARN "No terms were found for this aspect with a corrected P-value <= $cutoff.\n";
	} else {
	    my @lastresort;
	    for my $pv ( @pvalues ) {
		if( $pv->{CORRECTED_PVALUE} <= $cutoff) {
		    my @paths = $pv->{NODE}->pathsToRoot;
		    for my $path ( sort { scalar @$a <=> scalar @$b } @paths ) {
			if( @$path < 3 ) {
			    push @{$lastresort[@$path]}, [ $pv->{NODE}->term, sprintf("p=%.3f",$pv->{CORRECTED_PVALUE})];
			} else {
			    $go_seen{$termFinder->aspect}->{$pv->{NODE}->term} = sprintf("p=%.3f",$pv->{CORRECTED_PVALUE});
			}
			last;
		    }
		}
	    }
	    if( ! exists $go_seen{$termFinder->aspect} && @lastresort) {
		# just report the deepest node(s)
		$go_seen{$termFinder->aspect}->{$_->[0]} = $_->[1] for @{$lastresort[-1]};
	    }
	}
	print WARN "\n\n";
	print OUT join("\t",(map { my $x = $_;
				   join(";", 
					map { sprintf("%s:%s",
						      $_, $go_seen{$x}->{$_} ) }
					keys %{ $go_seen{$x} }) } qw(P C F) ),
		       ), "\n";
    }
    close(OUT);
    close(WARN);
}
