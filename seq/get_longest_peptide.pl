#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

# this finds the longest peptide when there are multiple isoforms 
# for a gene
# we assume gene names are encoded in the FASTA header
my $longest_ext = "longest";
my $min_length  = 50;

GetOptions(
	   'longest|ext:s'   => \$longest_ext,
	   'ml|min|minlen:i' => \$min_length,
	   );

for my $seqfile ( @ARGV ) {    
    my ($fh,$outfile);
    my ($stem,$ext,$cmp);
    if( $seqfile =~ /(\S+)\.(fa(?:s|sta)?)(?:\.(gz|bz2))?/ ) {
	($stem,$ext,$cmp) = ($1,$2,$3);
    } 
    if( defined $cmp ) {
	if( $cmp eq 'gz' ) {
	    open($fh,"zcat $seqfile |") || die $!;
	} elsif( $cmp eq 'bz2' ) {
	    open($fh,"bzip2 -c $seqfile |") || die $!;
	} else {
	    warn("unknown extension $ext!\n");
	}
    } else {
	open($fh, "< $seqfile") || die "cannot open $seqfile!";
    }
    $outfile = sprintf("%s.%s.%s",$stem,$longest_ext,$ext);
    
    my $in = Bio::SeqIO->new(-format => 'fasta',
			     -fh     => $fh);
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$outfile");
    my %seqs;
    my @order;
    while( my $seq = $in->next_seq ) {
	next if $seq->length < $min_length;
	if( $seq->description =~ /gene[=:](\S+)/ ) {
	    my $gn = $1;
	    push @order, $gn unless exists $seqs{$gn};
	    push @{$seqs{$gn}}, $seq;
	} else { 
	    warn("Cannot find a gene name encoded in fasta header: \n",
		 $seq->id, " ", $seq->description, "\n");
	    next;
	}
    }
    for my $gene ( @order ) {
	my ($longest) = sort { $b->length <=> $a->length } @{$seqs{$gene}};
	$out->write_seq($longest);
    }
}
