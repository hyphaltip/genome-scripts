#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SeqIO;
use Bio::Tools::Run::Alignment::TCoffee;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Align::Utilities qw(:all);
use Getopt::Long;
my $trans_only = 0;
my $aln_alg = 'muscle';
my $sp_sep = "_"; # was ":" in Stajich et al 2007
GetOptions(
	   't|trans'      => \$trans_only,
           'aln|alg:s'    => \$aln_alg,
	  );
my $force = 0; # whether or not to force regeneration of alignments
my ($factory,@params);
if( $aln_alg =~ /coffee/ ) { 
 @params = ('ktuple' => 2, 'matrix' => 'BLOSUM', -verbose => -1);
 $factory = Bio::Tools::Run::Alignment::TCoffee->new(@params);
} else {
 $factory = Bio::Tools::Run::Alignment::Muscle->new(@params);
}
$factory->quiet(1);
my $dir = shift || 'ortholog_sequence_sets';
my $outdir = "$dir-aln";
mkdir("$outdir") unless -d $outdir;
opendir(DIR, $dir) || die ($!);
my $testout = Bio::AlignIO->new(-format => 'clustalw');
FILE: for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.cds\.fa$/;
    
    my ($stem) = ($1);
    next if ( $trans_only && -e "$outdir/$stem.pep.fa");
    next if( ! $force && -e "$outdir/$stem.cds.phy" && 
	     ! -z "$outdir/$stem.cds.phy");
    my %cds_seqs;
    my $seqin = Bio::SeqIO->new(-format => 'fasta',
				-file   => "$dir/$file");
    my @prot_seqs;
    while( my $s = $seqin->next_seq ) {
	my $id = $s->display_id;
	$s->seq(uc($s->seq));
	my $pepseq = $s->translate();
	if( $pepseq->seq() =~ /\*$/ ) {
	    $pepseq = $pepseq->trunc(1,$pepseq->length-1);
	    $s = $s->trunc(1,$s->length - 3);
	    die $file unless( $s->translate->seq eq $pepseq->seq );	    
	}
	
	if( $pepseq->seq =~ /\*/ ) {
	    my $mod = ($s->length % 3);
	    if( $mod != 0 ) {
		$s = $s->trunc($mod+1,$s->length);
		$pepseq = $s->translate;
	    } else {
		warn("gene sequence '$id' in '$file' contains stop codons\n");
		next FILE;
	    }
	    if( $pepseq->seq =~ /\*$/ ) {
		$pepseq = $pepseq->trunc(1,$pepseq->length-1);
		$s = $s->trunc(1,$s->length - 3);
		die $file unless( $s->translate->seq eq $pepseq->seq );
	    } 
	    if( $pepseq->seq =~ /\*/ ) {
		warn("gene sequence '$id' in '$file' contains stop codons\n");
		next FILE;
	    }
	}

	if( $id =~ /^(\S+)[:\|]/) { $id = $1;}
	$s->display_id($id);
	$pepseq->display_id($id);
	
	$cds_seqs{$id} = $s;
	push @prot_seqs, $pepseq;
    }
    my $aln;    
    if( -e "$outdir/$stem.pep.fasaln" ) {
	my $in = Bio::AlignIO->new(-format => 'fasta',
				 -file   => File::Spec->catfile
				 ($outdir,"$stem.pep.fasaln"));
	unless( $aln = $in->next_aln ) {
	    warn("could not parse alignment for $stem.pep.fasaln\n");
	    next FILE;
	}
    } else {
	my $o = Bio::SeqIO->new(-format => 'fasta',
				-file   => '>'.File::Spec->catfile
				($outdir,"$stem.pep.fa"));
	$o->write_seq(@prot_seqs);
	next if $trans_only;
	$aln = $factory->align(\@prot_seqs);
    }

    $aln->set_displayname_flat(1); # simple output
    # fix aln
    {
	my @s = $aln->each_seq;
	for my $seq ( @s ) {
	    my $id = $seq->display_id;
	    if( $id =~ /^(\S+)[:\|]/) { $id = $1;}	    	
	    $aln->remove_seq($seq);
	    $seq->display_id($id);
	    $aln->add_seq($seq);
	}
    }
    my $cdsaln = &aa_to_dna_aln($aln,\%cds_seqs);
    if( ! $cdsaln->is_flush ) {
	$testout->write_aln($cdsaln);
	warn("trouble with $file\n");
	last;
    }
    my $outn = Bio::AlignIO->new(-format => 'nexus',
				-file   => ">$outdir/$stem.cds.nex");
    my $out = Bio::AlignIO->new(-format => 'phylip', -interleaved => 0,
				-idlength => 12,
				-file   => ">$outdir/$stem.cds.phy");
    $outn->write_aln($cdsaln);
    if( ! -e "$outdir/$stem.pep.fasaln" ||
	-z "$outdir/$stem.pep.fasaln" ) {
	my $outp = Bio::AlignIO->new(-format => 'fasta',
				     -file   => ">$outdir/$stem.pep.fasaln");
	$outp->write_aln($aln);    
    }
    $out->write_aln($cdsaln);
}
