#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::DB::Fasta;
use Getopt::Long;
use File::Spec;

my $pepdb = 'pep';
my $cdsdb = 'cds';
my $orth_file;

my $pepaln_dir = 'orthologs_pep';
my $cdsaln_dir = 'orthologs_cds_aln';
my $pepext     = 'pep.probcons.aln';
my $cdsext     = 'cds.aln';
my $aformat    = 'fasta';
my $oformat    = 'phylip';
my $prefix     = 'orthmcl';
my %extra_format = ( 'phylip' => [-interleaved => 1],
		     'nexus' => [-show_symbols => 0,
				 -show_endblock => 0]);

GetOptions(
	   'i|pd:s'          => \$pepaln_dir,
	   'o|cd:s'          => \$cdsaln_dir,
	   'pe:s'            => \$pepext,
	   'ce:s'            => \$cdsext,
	   'pepdb:s'         => \$pepdb,
	   'cdsdb:s'         => \$cdsdb,
	   'aformat:s'       => \$aformat,
	   'oformat:s'       => \$oformat,	   
	   'prefix:s'        => \$prefix,
	   );
mkdir($cdsaln_dir) unless -d $cdsaln_dir;
my $cds_db = Bio::DB::Fasta->new($cdsdb);

while(<>) {
    my ($id,@genes) = split;
    my $pepalnfile = File::Spec->catfile($pepaln_dir, "$prefix\_$id.$pepext");
    my $cdsalnfile = File::Spec->catfile($cdsaln_dir, "$prefix\_$id.$cdsext");
    if( ! -f $pepalnfile ) {
	warn("no alignment for $pepalnfile  in $pepaln_dir\n");
    
	next;
    }
    my $in = Bio::AlignIO->new(-format => $aformat,
			       -file   => $pepalnfile);
    my $pepaln = $in->next_aln;
    if( ! $pepaln ) {
	warn("no alignment found in file $pepalnfile\n");
	next;	
    }    
    my %seqs;
    for my $gene (@genes ) {
	my ($id) = split(/:/,$gene);

	my $geneseq = $cds_db->seq($gene);
	if( ! $geneseq ) {
	    warn("no sequence for $gene found in db $cdsdb\n");
	}
	$seqs{$id} = Bio::Seq->new(-seq => uc $geneseq, 
				   -id  => $gene);
    }
    my $dna_aln = &aa_to_dna_aln($pepaln,\%seqs);
    my $out = Bio::AlignIO->new(-format => $oformat,
				-file   => ">$cdsalnfile",
				@{$extra_format{$oformat} || []}
				);
    for my $seq ( $dna_aln->each_seq ) {
      my $id = $seq->display_id;
      $id =~ s/[\.\-]/\_/g;
      $dna_aln->remove_seq($seq);
      $seq->display_id($id);
      $dna_aln->add_seq($seq);
    }
    $out->write_aln($dna_aln);
}
