#!/usr/bin/perl -w
use strict;
# we want this
# scaffold_1 src exon 9975 11074   .  +   . gene_id S.10001 transcript_id S.10001 exon_number 1
# scaffold_1  src  CDS 10027 10824   .  +   . gene_id S.10001 transcript_id S.10001 exon_number 1
my %genes;
my @mrna;
while(<>) {
    next if /^\#/;
    chomp;
    my @line = split(/\t/,$_);
    my $last_col = $line[-1];
    my %group = map { split(/=/,$_) } split(/;/,pop @line);
    my $group = $group{'Parent'};
    if( $line[2] eq 'CDS' || $line[2] eq 'exon' ) {
        push @mrna, $group unless exists $genes{$group};
	push @{$genes{$group}}, [@line];
    }
}

for my $gene ( @mrna ) { # by gene here, I really mean transcript_name
    my $exonct = 1;
    my @x;
    if( ! exists $genes{$gene} ) {
     warn("no value for $gene\n");
     next;
    }
    my $genename= $gene;
   $genename =~ s/\.t\d+$//;
	
    my @exons = sort { $a->[3] * ($a->[6] eq '-' ? -1 : 1) <=>
                              $b->[3] * ($b->[6] eq '-' ? -1 : 1) }
                   @{$genes{$gene}} ;
    for my $exon ( @exons ) {
	push @x, [@$exon, sprintf("gene_id %s transcript_id %s exon_number %d",
			      $genename, $gene, $exonct)];
	$exonct++ if $exon->[2] eq 'CDS';
        
    }
    for my $exon ( sort { $a->[3] <=> $b->[3] } @x ) {
	print join("\t", @$exon), "\n";
    }
}
