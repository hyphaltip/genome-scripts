#!/usr/bin/perl -w
use strict;

my $coprinus_genes = 'FISH/coprinus.simple.gff3';
my $laccaria_genes = 'FISH/laccaria.simple.gff3';
my $fish = 'FISH/coprinus_laccaria_FISH.txt';

my %num2roman = ( '1' => 'I',
		  '2' => 'II',
		  '3' => 'III',
		  '4' => 'IV',
		  '5' => 'V',
		  '6' => 'VI',
		  '7' => 'VII',
		  '8' => 'VIII',
		  '9' => 'IX',
		  '10'=> 'X',
		  '11'=> 'XI',
		  '12'=> 'XII',
		  '13'=> 'XIII',
		  );
my %gene2pos;
open(my $fh => "cat $coprinus_genes $laccaria_genes |") || die $!;
while(<$fh> ) {
    next if /^\#/;
    chomp;
    my @row = split(/\t/,$_);
    my %last = map { split(/=/,$_) } split(/;/,pop @row);
    $gene2pos{$last{ID}} = [ $row[0], $row[3],$row[4] ];
}
close($fh);

open($fh, $fish) || die $!;
my @blocks;
my $block = 0;
while(<$fh>) {
    if( /^(\d+)\s+$/ ) {
	$block = $1;
    } elsif( /^\s+$/) {
	next;
    } else {
	my @row = split;
# col1:   Genome1
# col2:   Chromosome(Scaffold) of Genome1
# Col3:   Block start gene index on the Chromosome or Scaffold of Genome1
# Col4:   Block end gene index of Genome1
# Col5:   Block start gene ID of Genome1
# Col6:   Block end gene ID of Genome1

# Col7:   Genome2
# Col8:   Chromosome(Scaffold) of Genome2
# Col9:   Block start gene index on the Chromosome or Scaffold of Genome2
# Col10:  Block end gene index of Genome2
# Col11:  Block start gene ID of Genome2
# Col12:  Block end gene ID of Genome2

# Col13:  Anchor gene ID from Genome1 (anchor gene is the homologous 
#          gene pair anchoring two syntenic blocks)
# Col14:  Anchor gene ID from Genome2
# Col15:  Anchor gene index from Genome1
# Col16:  Anchor gene index from Genome2
	next unless $row[0] ne $row[6];
	#print (join("\t", $block,$row[4], $row[5], $row[10], $row[11]),"\n");
	$blocks[$block] = [ $row[4], $row[5], $row[10], $row[11] ];
    }
}
my $i = 0;
print '#',join("\t", qw(chrom start stop block)),"\n";
for my $block ( @blocks ) {
    if( defined $block ) {
	my ($gene_start) = $gene2pos{$block->[0]};
	my ($gene_stop) = $gene2pos{$block->[1]};
	my $chr = $gene_start->[0];

	$chr =~ s/Chrom(\d+)/$num2roman{$1}/;
	print join("\t", $chr,  # chrom
		   $gene_start->[1],        # block start, gene start
		   $gene_stop->[2],         # block stop, gene stop
		   sprintf("Block_%d",$i)),"\n";
    }
    $i++;
}

