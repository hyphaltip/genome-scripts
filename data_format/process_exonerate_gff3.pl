#!/usr/bin/perl -w

=head1 NAME 

process_exonerate_gff3 -t EST file1.e2g.EXONERATE file2.e2g.EXONERATE > out.gff3

=head1 USAGE

Command line arguments
 -gf|gff_version  gff version string (1,2,2.5,3) default 3
 -t|type EST or Protein typically, the string in the ID= field default EST 

=head1 DESCRIPTION

Turns EXONERATE gff output into GFF for Gbrowse use.

You need to have run exonerate with at least the following options

 --showtargetgff yes
 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"

Although I often use this for fungal protein mapping

 --showvulgar yes --softmaskquery yes --softmasktarget yes --minintron 20 --maxintron 3000 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --showalignment no

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use IO::String;
use Getopt::Long;
use Bio::Tools::GFF;
use Bio::Seq;
use Env;
use File::Spec;
my $state = 0; my $ct = 0;
my $buffer = '';
my (%dat,%lengths,%alnlens,%counter,$vulgar,$lastvulgar);

my $type = 'EST'; # or Protein
my $gff_ver = 3;
my $show_aln = 0;
GetOptions(
	   'a|showaln!'          => \$show_aln,
	   't|type:s'            => \$type,
	   'gff_version|gf|v:s' => \$gff_ver,
           'h|help'              => sub { exec('perldoc', $0); exit }
	   );
my $out = Bio::Tools::GFF->new(-gff_version => $gff_ver);


while(<>) {
    if( $state == 0 && s/^vulgar:\s+// ) {
	chomp($vulgar = $lastvulgar = $_);
    }
    if( $state > 0 ) {	 
	if( $state == 2 ) {
	    if ( s/^vulgar:\s+//) {
		$lastvulgar = $vulgar;
		chomp($vulgar = $_);
	        chomp($lastvulgar);
	    }
	    my $in = Bio::Tools::GFF->new
		(-gff_version => 2,
		 -fh          => IO::String->new($buffer));
	    my ($gene_name,$gene);
	    my $length;
	    while( my $f = $in->next_feature ) {
	        my $srctag  = $f->source_tag;
		$srctag =~ s/\:/_/g;
		$f->source_tag($srctag);
		if( $f->primary_tag eq 'gene' ) {
		    $length = $lengths{$f->seq_id};		    
		    if(  ! defined $length ) {
			die("unknown length for '",$f->seq_id,"'\n");
		    }
		    ($gene) = $f->get_tag_values('sequence');
		    if ( $f->has_tag('gene_orientation')) {
			if( $dat{$gene} ) {
			    $gene .= $counter{$gene};
			}
			($dat{$gene}) = $f->get_tag_values('gene_orientation');
			$dat{$gene} = 0 if $dat{$gene} eq '.';
			$f->remove_tag('gene_orientation'); 
		    }
		    for my $t ( qw(gene_id sequence) ) {
			$f->remove_tag($t) if $f->has_tag($t); 
		    }
		    $gene_name = $gene;
		    if( $counter{$gene_name}++ ) {
			$gene_name .= ".$counter{$gene_name}";
		    }

		    $f->add_tag_value('ID',"$type:$gene_name");
		    chomp($lastvulgar);
		    $f->add_tag_value('vulgar', $lastvulgar);
		} elsif( $f->primary_tag eq 'similarity') {
		    next unless $show_aln;
		    if( $f->has_tag('Align') ) { 
			my @hsps = $f->get_tag_values('Align');
			$f->primary_tag('match');
			$f->remove_tag('alignment_id');
			
			my ($min,$max);
			while ( @hsps ) {
			    my ($hstart, $qstart,$hlen) = splice(@hsps,0,3);
			    my $hlen_mod = $hlen;
			    if( $type =~ /Protein/i ) {
				$hlen_mod /= 3;
			    }
			    my $newf = Bio::SeqFeature::Generic->new
				( -seq_id => $f->seq_id,
				  -start => $hstart,
				  -end   => $hstart + $hlen,
				  -strand=> $dat{$gene},
				  -source_tag => $srctag,
				  -primary_tag => 'HSP',
				  -tag => {
				      'Target' => sprintf("$type:%s %d %d",
							  $gene_name,
							  $qstart,
							  $qstart+$hlen_mod),
				  }
				  );
			    $out->write_feature($newf);
			    $min = $qstart unless defined $min;
			    $max = $qstart + $hlen_mod;

			}
			$f->add_tag_value('Target',
				    sprintf("$type:%s %d %d",
					    $gene_name, $min, $max,
					    ));
			my $genelen = $lengths{$gene};
			$f->add_tag_value('Length', $genelen);
			$f->add_tag_value('Coverage', int(100 * ($max-$min)/$genelen));  
			$f->remove_tag('Align');
		    }
		} elsif( $f->primary_tag eq 'splice5' ) {
		    $f->primary_tag('splice_donor');
		} elsif( $f->primary_tag eq 'splice3' ) {
		    $f->primary_tag('splice_acceptor');
		} elsif( $f->primary_tag eq 'frameshift' ) {
		    $f->primary_tag('frameshift');
		} elsif( $f->primary_tag eq 'exon' ) {
		    $f->primary_tag('CDS');
		}
		for my $tag ( qw(intron_id deletions insertions) ) {
		    $f->remove_tag($tag) if $f->has_tag($tag); 
		}
		if( $f->strand < 0 ) {
		    my $s = $length - $f->end + 1;
		    my $e = $length - $f->start +1;
		    $f->start($s);
		    $f->end  ($e);
		}
		if( $gene && 
		    $dat{$gene} &&
		    $dat{$gene} eq '-' ) 
		{
		    $f->strand(-1);
		}
		if( $f->primary_tag ne 'gene' && $f->primary_tag ne 'match') {
		    $f->add_tag_value('Parent', "$type:$gene_name");
		}
		$out->write_feature($f);
	    }
	    $state = 0;
	    $buffer = '';
	    %dat = ();
	    
	    next;
	} elsif(/^\#/) { next; }
	if( /^>(\S+)\s+length=(\d+)\s+alnlen=(\d+)/ ) {
	    $lengths{$1} = $2;
	    $alnlens{$1} = $3;
	    $ct++;
	} else {	
	    $buffer .= join("\t",split(/\s+/,$_,9));
	}
	if( $ct == 2 ) { 
	    $state = 2;
	    $ct = 0;
	}
    } elsif( /^\# --- START OF GFF DUMP ---/ ) {
	$state = 1; 
	$ct =0;
	$buffer = '';
    }
}
