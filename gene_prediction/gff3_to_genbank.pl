# -*-perl-*-
use strict;
use warnings;

use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::DB::SeqFeature;

my $gff = shift;
my $fasta = shift;

my $out = Bio::SeqIO->new(-format => 'genbank');

my $in = Bio::SeqIO->new(-format => 'fasta', -file => $fasta);
while(my $f = $in->next_seq ) {
    $seqs{$f->display_id} = $f->seq;
}

my %seqs;
my %feats;
open( my $fh => $gff) || die $!;
my %cds2gene;
my %mrna2gene;
while(<$fh>) {
    next if /^\#/;
    chomp;
    my @row = split(/\t/,$_);
    my $id = $row[0];
    next if $row[2] !~ /gene|mRNA|CDS/;
    my %group = map { split(/=/,$_) } split(/;/,pop @row);
    my $group_id = $row[2] eq 'CDS' ? $group{'Parent'} : $group{'ID'};
    my ($group_id,$gene_id);
    if( $row[2] eq 'CDS' ) {
	$group_id = $group{'Parent'};
	$gene_id = $mrna2gene{$group_id};
    } elsif( $row[2] eq 'mRNA' ) {
	$group_id = $group{'ID'};
	$gene_id = $group{'Parent'};
    } elsif( $row[2] eq 'gene' ) {
	$gene_id = $group{ID};
    } else {
	warn("skilling $row[2] ",$group{ID},"\n");
    }
	
    push @{$feats{$id}->{$row[2]}->{$group_id}
	
}
for my $ctg ( keys %feats ) {
    my $seq = $seqs{$ctg};
    for my $g ( @{$ctgs{$ctg}} ) {
	
	my $gene = Bio::SeqFeature::Generic->new(-start => $g->
	$seq->add_SeqFeature($gene);
    }
}
