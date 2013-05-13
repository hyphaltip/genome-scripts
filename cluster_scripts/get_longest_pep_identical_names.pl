#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

=head1 NAME

get_longest_by_name_identical - get longest sequence when multiple exist with same name

=head1 SYNOPSIS

perl get_longest_by_name_identical pepfile > pepfile.nr

=head1 DESCRIPTION

Take the longest sequence when there are multiple with same name

=head1 AUTHOR 

Jason Stajich, jason.stajich[at]ucr.edu

=cut


my %seqs;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => shift);
my $out = Bio::SeqIO->new(-format => 'fasta');


while( my $s = $in->next_seq ) {
 if( ! exists $seqs{$s->display_id} ) {
   $seqs{$s->display_id} = $s;
 } else {
  if( $seqs{$s->display_id}->length < $s->length) {
   $seqs{$s->display_id} = $s;
  }
 }
}
for my $nm ( sort keys %seqs ) {
 $out->write_seq($seqs{$nm});
}
