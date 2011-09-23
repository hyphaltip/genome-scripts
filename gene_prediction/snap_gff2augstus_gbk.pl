#!/usr/bin/perl -w
use strict;
use List::Util qw(min max);
use Getopt::Long;
use Bio::SeqIO;
#use Bio::DB::GFF;
use Bio::SeqFeature::Generic;
use Bio::Location::Split;
use Bio::Location::Simple;
use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;
use File::Temp qw(tempdir tempfile);

my ($gff,$fa,$out,$target);

GetOptions(
	   'i|in|gff:s'  => \$gff,
	   'fa|s|seq:s'  => \$fa,
	   'o|out:s'     => \$out,
	   't|target:s'  => \$target,
	   );

if( ! defined $gff ) {
    die "must provide a gff file with 'i|in|gff'\n";
}
if( ! defined $fa ) {
    die "must provide a fasta file with 'fa|s|seq'\n";
}

my $out_seq;
if( defined $out ) {
    $out_seq = Bio::SeqIO->new(-format => 'genbank',
			       -file   => ">$out");
} else {
    $out_seq = Bio::SeqIO->new(-format => 'genbank');
}
my $tmpdir = tempdir(CLEANUP => 1);
#my $db = Bio::DB::GFF->new(-adaptor => 'berkeleydb',
#			   -create  => 1,
#			   -dsn     => $tmpdir,
#			    );
my ($fh,$gff_chroms,$db);
if( $target ) {
    ($gff_chroms) = "$target.gff";

    if( -d $target ) {
	$db = Bio::DB::SeqFeature::Store->new(-adaptor=> 'berkeleydb3',
					      -fa     => $fa,
					      -dsn    => $target);
    } else {
	$db = Bio::DB::SeqFeature::Store->new(-adaptor=> 'berkeleydb3',
					      -dsn    => $target,
					      -fasta  => $fa,
					      -create => 1);
	open($fh => ">$gff_chroms") || die $!;
    }
} else {
    $db = Bio::DB::SeqFeature::Store->new(-adaptor=> 'berkeleydb3',
					  -dsn    => $tmpdir,
					  -fasta  => $fa,
					  -create => 1);
    ($fh,$gff_chroms) = tempfile(); 	
}
if( $fh ) {
    print $fh "##gff-version 3\n";
    my @seqs;
    my $infa = Bio::SeqIO->new(-format => 'fasta', -file => $fa);

    while(my $seq = $infa->next_seq ) {
	my $nm = $seq->display_id;
	print $fh join("\t", $nm, 'chromosome', 'scaffold',
		       1, $seq->length, '.', '.', '.', 
		       sprintf("ID=%s;Name=%s",$nm,$nm)
		       ), "\n";
	push @seqs, $seq;
    }
    
    open(my $in => $gff ) || die $!;
    my %genes;
    my (%len);
    my $last_nm;
    while(<$in>) {
	next if /^\#/;
	chomp;
	my (@row) = split(/\t/,$_);
	if( defined $last_nm && $last_nm ne $row[0] ) {
	    print $fh join("\t", $last_nm, qw(snap mRNA),
			   $len{$last_nm}->[0],
			   $len{$last_nm}->[1],'.',
			   '+', '.',
			   sprintf("ID=%s_mRNA;Parent=%s_gene",$last_nm,$last_nm)),"\n";
	    print $fh join("\t", $last_nm, qw(snap gene),
			   $len{$last_nm}->[0],
			   $len{$last_nm}->[1],'.',
			   '+', '.',
			   sprintf("ID=%s_gene",$last_nm)),"\n";
	}
	$row[-1]  = '.'; # missing the frame column;
	push @row, sprintf("Parent=%s_mRNA",$row[0],$row[0]);
	print $fh join("\t", @row), "\n";
	
	$len{$row[0]}->[0] = defined $len{$row[0]}->[0] ?  
	    min($row[3], $len{$row[0]}->[0]) :
	    $row[3];
	$len{$row[0]}->[1] = defined $len{$row[0]}->[1] ?  
	    max($row[4], $len{$row[0]}->[1]) :
	    $row[4];
	$last_nm = $row[0];
    }
    
    if( $last_nm ) {
	print $fh join("\t", $last_nm, qw(snap mRNA),
		       $len{$last_nm}->[0],
		       $len{$last_nm}->[1],'.',
		       '+', '.',
		       sprintf("ID=%s_mRNA;Parent=%s_gene",$last_nm,$last_nm)),"\n";
	print $fh join("\t", $last_nm, qw(snap gene),
		       $len{$last_nm}->[0],
		       $len{$last_nm}->[1],'.',
		       '+', '.',
		       sprintf("ID=%s_gene",$last_nm)),"\n";
	
	
    }
    
    close($fh);
# Load an entire GFF3 file, using the GFF3 loader...
    my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $db,
							     -verbose  => 1,
							     -fast     => 1);
    warn($gff_chroms, "\n");
    $loader->load($gff_chroms);
    $loader->load($fa);
    unlink($gff_chroms) unless $target;
}

my $iterator = $db->get_seq_stream(-type => 'gene');

while( my $gene = $iterator->next_seq ) {
    my $chrom_seq = Bio::Seq::RichSeq->new(-seq => $db->seq($gene->seq_id),
					   -accession_number => $gene->seq_id,
					   -display_id => $gene->seq_id,
					   );
    my $src_feature = Bio::SeqFeature::Generic->new
	(-primary_tag => 'source',
	 -start       => 1,
	 -end         => $chrom_seq->length);
    $chrom_seq->add_SeqFeature($src_feature);

    for my $mRNA ( $gene->get_SeqFeatures('mRNA') ) {
	my $loc;
	my @cds = $mRNA->get_SeqFeatures('CDS');
	if( @cds == 1 ) {
	    $loc = $cds[0]->location;
	} else {
	    $loc = Bio::Location::Split->new();
	    for my $cds ( sort { $a->start <=> $b->start } @cds ) {
		$cds->location->seq_id('');
		$loc->add_sub_Location($cds->location);
	    }
	}
	my $cds_feat = Bio::SeqFeature::Generic->new
	    (-primary_tag => 'CDS',
	     -location    => $loc);
	$chrom_seq->add_SeqFeature($cds_feat);
						     
    }
    $out_seq->write_seq($chrom_seq);
}
