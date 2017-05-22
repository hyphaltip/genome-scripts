#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::AlignIO; # bioperl

=head1 NAME

aln_identify_unique_cols.pl - identify unique residue/nt configurations

=head1 DESCRIPTION

For a given alignment print out classification for each column and the 
uniqueness of the base for each column, given a target sequence prefix/id

=head1 USAGE

Typical usage is where sequences in alignments are prefixed with a
species name, e.g. Scer|YAL001 - the prefix provided would be Scer


perl aln_identify_unique_cols.pl -p PREFIX -f fasta alignment.aln

=head1 AUTHOR

Jason Stajich jasonstajich.phd[at]gmail.com

=cut

my $prefix;
my $format = 'fasta';
GetOptions('ref|p|prefix:s' => \$prefix,
	   'f|format:s'     => \$format,
    );


if ( ! defined $prefix ) {
    die("requires a target prefix for the reference to consider\n");
}
my $alnfile = shift || die "need align file";

my $alnio = Bio::AlignIO->new(-file   => $alnfile,
			      -format => $format);

while(my $aln = $alnio->next_aln ) {
    $aln->set_displayname_flat(1);
    my @matrix;
    my $ref_pos = -1;
    my $seqn = 0;
    for my $seq ( $aln->each_seq ) {
	my ($pref,$id) = split(/\|/,$seq->display_id);
	if( $pref eq $prefix ) {
	    warn("found target seq $pref -> $id\n");
	    $ref_pos = $seqn;
	}
	my $col = 0;
	for my $c ( split('',$seq->seq) ) { 
	    $matrix[$col++] .= $c;
	}
	$seqn++;
    }
    my $i =0 ;
    
    for my $row ( @matrix ) {

	my %let_freq;
	for my $bp ( split('',$row) ) {
	    $let_freq{$bp}++;
	}
	my $info = '';
	my @freq = keys %let_freq;
	if( scalar @freq == 1 ) {
	    $info = 'Invariant';
	} elsif ( scalar @freq == 2 ) {
	    my $ref_base = substr($row,$ref_pos,1);
	    if( $let_freq{$ref_base} == 1 ) {
		# only 2 bases, and ref base is different
		$info = sprintf("Target (%s) is Uniquely Different",$prefix);
	    }
	}
		
	printf "%-4d => %s %s\n",$i++,$row, $info;	   
    }
}
