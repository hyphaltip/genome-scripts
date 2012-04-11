#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Bio::DB::Fasta;

# Currently this assume COMPLETE PASA calls

my ($keep, $db, $gff,$prefix);

my $suffix = 'train';

GetOptions(
	   'i|in|keep:s' => \$keep, # list of assemblies to use
	   'g|gff:s'     => \$gff, # gff PASA result
	   'd|db:s'      => \$db, # genomic db
	   'p|prefix:s'  => \$prefix,
	   );
die("need a list of transcripts to use with -i or --keep") unless defined $keep;
die("need a PASA gff file with -g or --gff") unless defined $gff;
die("need a db with -d or --db") unless defined $db;
die("need a prefix with -p or --prefix") unless defined $prefix;
my $dbh = Bio::DB::Fasta->new($db);

my %targets;

open( my $fh => $keep) || die "$keep: $!";
while(<$fh>) {
    my ($name) = split;
    $name =~ s/^>//;
    $targets{$name}++;
}
close($fh);
open(my $gffh => $gff) || die $!;
# PASA generates GFF3

my %keep;
while(<$gffh>) {
    next if /^\s+$/;
    chomp;
    my @line = split(/\t/,$_);
    next unless $line[2] eq 'CDS';
    my %group = map { split(/=/,$_,2) } split(/;/,pop @line);
    my $parent = $group{'Parent'};
    $parent =~ s/LongOrf\.//;
    if( exists $targets{$parent} ) {
	my ($s,$e) = ($line[3],$line[4]);
	if( $line[6] eq '-' ) {
	    ($e,$s) = ($s,$e);
	}
	push @{$keep{$line[0]}->{$parent}}, [$s,$e,$parent];
    }
}
close($gffh);
my $outseq = Bio::SeqIO->new(-format => "fasta",
			     -file   => ">$prefix.train.fa");
open($fh => ">$prefix.train.zff") || die $!;
while( my ($chrom,$genes) = each %keep ) {
    printf $fh ">%s\n",$chrom;
    for my $gene ( keys %{$genes} ) {
	my $exons = $genes->{$gene};
	my $ecount = scalar @$exons;
	my $i = 1;
	for my $exon ( @$exons ) {
	    my $exontype = 'Exon';
	    if( $i == 1 ) {
		if( $i == $ecount ) {
		    $exontype = 'Esngl';
		} else {
		    $exontype = 'Einit';		    
		}
	    } elsif( $i == $ecount ) {
		$exontype = 'Eterm';
	    }
	    print $fh join("\t",$exontype, @{$exon}),"\n";
	    $i++;
	}
    }
    $outseq->write_seq($dbh->get_Seq_by_acc($chrom));
}
