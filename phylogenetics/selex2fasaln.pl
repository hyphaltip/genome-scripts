#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;

# used to convert hmmalign results to input to phylogenetic applications

my $iformat = 'selex';
my $oformat = 'fasta';
my $file = shift;
my $in = Bio::AlignIO->new(-format => $iformat,
			   -file     => $file);
my $out = Bio::AlignIO->new(-format => $oformat);
while( my $aln = $in->next_aln ) {
    $aln->map_chars('[\*\.]','-');
    my $grm = $aln->remove_gaps(undef,1);
    $grm->set_displayname_flat(1);
    $out->write_aln($grm);
}
