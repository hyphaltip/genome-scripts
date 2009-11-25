#!/usr/bin/perl -w
use strict;
use IO::String;
use Bio::Tools::GFF;
use Bio::DB::Fasta;
use Env;
use File::Spec;
my $out = Bio::Tools::GFF->new(-gff_version => 3);

my $db = Bio::DB::Fasta->new(File::Spec->catfile($HOME,qw(fungi 
							  fungal_genomes
							  nt)));
my $state = 0;
my $buffer = '';
while(<>) {

    if( $state == 1 ) {
	if( /^\# --- END OF GFF DUMP ---/) { 	    
	    my $in = Bio::Tools::GFF->new(-gff_version => 2,
					  -fh          => 
					  IO::String->new($buffer));
	    my $gene;
	    my $length;
	    while( my $f = $in->next_feature ) {
	        my $srctag  = $f->source_tag;
		$srctag =~ s/\:/_/g;
		$f->source_tag($srctag);
		if( $f->primary_tag eq 'gene' ) {
		    $length ||= $db->length($f->seq_id);
		    ($gene) = $f->get_tag_values('sequence');
		    $f->remove_tag('gene_orientation') if $f->has_tag('gene_orientation');
		    $f->add_tag_value('Id',$gene);
		    if( $f->strand < 0 ) {
			# we're flip-flopping start/end
			my $s = $length - $f->end;
			my $e = $length - $f->start;
			$f->start($s);
			$f->end  ($e);
		    }
		    $out->write_feature($f);
		    next;
		} elsif( $f->primary_tag eq 'similarity') {
		    # make sub pieces
		    next;
		    #$f->add_tag_value('Target', $gene);
		} elsif( $f->primary_tag eq 'splice5' ) {
		    $f->primary_tag('splice_donor');
		} elsif( $f->primary_tag eq 'splice3' ) {
		    $f->primary_tag('splice_acceptor');
		} elsif( $f->primary_tag eq 'exon' ) {
		    $f->primary_tag('CDS');
		}
		if( $f->strand < 0 ) {
		    my $s = $length - $f->end;
		    my $e = $length - $f->start;
		    $f->start($s);
		    $f->end  ($e);
		}
		$f->add_tag_value('Parent', $gene);
		$out->write_feature($f);
	    }
	    $state = 0;
	    $buffer = '';

	    next;
	} elsif(/^\#/) { next; }
	$buffer .= join("\t",split(/\s+/,$_,9));
	
    } elsif( /^\# --- END OF GFF DUMP ---/) { 
	$state = 0;
	$buffer = '';
    } elsif( /^\# --- START OF GFF DUMP ---/ ) {
	$state = 1;
	$buffer = '';
    } 
}
