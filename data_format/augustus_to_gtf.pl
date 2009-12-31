#!/usr/bin/perl -w
use strict;

# Convert AUGUSTUS GFF3 to GTF suitable for
# input to most gene prediction programs for training
# or (for our purposes) as input to Mercator

use Getopt::Long;
my $prefix = 'augustus';

GetOptions(
	   'p|prefix=s' => \$prefix,
	   );
my (%codons,@exons);
my $gene;
while(<>) {
    if( /\# gene (\S+)/ ) {
	$gene = $1;	
    } elsif(/\# end gene /) {
	if( exists $codons{'start_codon'} ) {
	    print join("\t",@{$codons{'start_codon'}}),"\n";
	}
	my $i = 0;
	my $exontype = (@exons == 1) ? 'single' : 'initial';
	for my $l ( map { $_->[0] }
		    sort { $a->[1] <=> $b->[1] }
		    map { [$_, $_->[3] * ($_->[6] eq '+' ? 1 : -1)] }
		    @exons ) {
	    if( $i == 0 && $l->[7] != 0 ) {
		if( $l->[6] eq '+') {
		    $l->[3] += $l->[7];
		} else {
		    $l->[4] -= $l->[7];
		}
	    }
	    $l->[-1] = sprintf("exontype \"%s\"; %s", $exontype,$l->[-1]);
	    print join("\t", @$l),"\n";
	    $l->[2] = 'exon';
	    print join("\t", @$l),"\n";
	    $exontype = ++$i == $#exons ? 'terminal' : 'internal';
	}
	if( exists $codons{'stop_codon'} ) {	
	    print join("\t",@{$codons{'stop_codon'}}),"\n" if exists $codons{'stop_codon'};
	}
	%codons = ();
	@exons  = ();
    }
    next if /^Make|Fehler/ || /^\s+$/ || /^\#/ || /addGene/;
    chomp;
    my @line = split(/\t/,$_);
    my $lastcol = pop @line;
    if( ! defined $line[2] ) {
	warn("no line[2] for $_");
	next;
    } elsif( $line[2] eq 'CDS' ) {
	$lastcol =~ s/ID=\S+;Parent=((\S+)\.t\d+)/transcript_id "$1"; gene_id "$2";/;
	push @exons, [@line, $lastcol];
    } elsif( $line[2] eq 'transcript' ) {	
	$line[2] = 'mRNA';
	
#	$lastcol =~ s/ID=(\S+);\S+/GenePrediction $line[0]-$prefix-$1/;
    } elsif( $line[2] eq 'gene' ) {
	next;
    } elsif( $line[2] eq 'start_codon' || 
	     $line[2] eq 'stop_codon' ) {
	$lastcol =~ s/Parent=((\S+)\.t\d+)/transcript_id "$2"; gene_id "$2";/;
	$codons{$line[2]} = [@line,$lastcol];
    } else { 
	next;
    }
}
