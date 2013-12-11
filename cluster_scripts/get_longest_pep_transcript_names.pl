#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

=head1 NAME

get_longest_pep_transcript_names - get longest sequence when there are multiple for same gene (alt splicing). 

=head1 SYNOPSIS

perl get_longest_pep_by_trasncript_name > pepfile.nr

=head1 DESCRIPTION

Take the longest sequence when there are multiple with same gene name, will attempt to parse gene name from the transcript name - so this mainly only works for T\d+ or \-R[A-Z] or \-R\d+ patterns

=head1 AUTHOR 

Jason Stajich, jason[at]bioperl.org

=cut


my %seqs;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => shift);
my $out = Bio::SeqIO->new(-format => 'fasta');


while( my $s = $in->next_seq ) {
    my $id = $s->display_id;
    my $gname = $id;
    if( $id =~ /(\S+)(T\d+|\-R([A-Z]|\d+))$/){
	$gname = $1;
    } else { 
	warn("cannot determine name for $id\n");
	next;
    }
    if( ! exists $seqs{$gname} ) {
	$seqs{$gname} = $s;
    } else {
	if( $seqs{$gname}->length < $s->length) {
	    $seqs{$gname} = $s;
	}
    }
}

for my $nm ( sort keys %seqs ) {
    $out->write_seq($seqs{$nm});
}
