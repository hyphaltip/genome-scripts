#!/usr/bin/perl -w
use strict;

use File::Spec;
use Getopt::Long;

use Bio::DB::SeqFeature;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::Location::Split;

my $dir;
my $ext = "gff3";
my $cdsdir = 'seqs/cds';
my $pepdir = 'seqs/pep';
my $introndir = 'seqs/intron';
my $genedir = 'seqs/gene';
my ($verbose,$force) = (0,0);
my $debug = 0;
GetOptions(
	   'f|force!'   => \$force,
	   'd|dir:s'    => \$dir,
	   'e|ext:s'    => \$ext,
	   'c|cds:s'    => \$cdsdir,
	   'p|pep:s'    => \$pepdir,
	   'i|intron:s' => \$introndir,
	   'g|gene:s'   => \$genedir,
	   'v|verbose'  => \$verbose,
	   'debug!'     => \$debug,
	   );
$dir = shift @ARGV unless $dir;
mkdir('seqs') unless -d 'seqs';
for my $d ( $cdsdir, $pepdir, $introndir, $genedir) {
    mkdir($d) unless -d $d;
}
die("need a dir with -d or as ARGV\n") unless defined $dir;
for my $dx ( $cdsdir, $pepdir, $introndir, $genedir ) {
    mkdir($dx) unless -d $dx;
}

opendir(DIR, $dir) || die $!;
for my $genome ( readdir(DIR) ) {
    next unless ( index($genome,".") != 0 && -d "$dir/$genome");
    
    my ($genus,$species,$strain,$version,$other) = split(/_/,$genome);
    if( $other ) {
	$strain = $version;	
    }
    if( ! $version ) {
	$version = $strain;
	$strain = undef;
    }
    my $prefix = substr($genus, 0,1) . substr($species,0,3);
    if( $strain ) {
	$prefix .= "_$strain";
    }
    
    warn("$genome $prefix\n");
    if( -f "$pepdir/$genome.pep.fa"  && ! -z "$pepdir/$genome.pep.fa" ) {
	unless( $force ) { 
	    warn("\tskipping $genome, already processed\n") if $verbose;
	    next;
	}
    }

    my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'berkeleydb',
#					     -dsn     => $dbdir,
					     -dir     => "$dir/$genome");

    my %types = $db->types;
    for my $t ( grep { /gene:/ } keys %types ) {
	next unless $t =~ /gene:(\S+)/;
	next if $t =~ /genewise/;
	my $base = $1;
	my $pout = Bio::SeqIO->new(-format => 'fasta',
				   -file   => ">$pepdir/$genome.$base.pep.fa");
    
	my $cout = Bio::SeqIO->new(-format => 'fasta',
				   -file   => ">$cdsdir/$genome.$base.cds.fa");
    
	my $iout = Bio::SeqIO->new(-format => 'fasta',
				   -file   => ">$introndir/$genome.$base.intron.fa");
	my $gout = Bio::SeqIO->new(-format => 'fasta',
				   -file   => ">$genedir/$genome.$base.gene.fa");
	
	my @features = $db->get_features_by_type($t);
	warn("there are ", scalar @features, " gene features\n") if $verbose;
	for my $f ( @features ) {
	    my ($gname) = $f->name || $f->load_id;
	    $gout->write_seq(Bio::PrimarySeq->new(-id => "$gname",
						  -seq=> $f->dna,
						  -desc => 
						  sprintf("%s:%s",$f->seq_id,
							  $f->location->to_FTstring())));
						  
	    my @mrna = $f->get_SeqFeatures('mRNA');
	    ### warn("there are ", scalar @mrna, " features for gene ", 
	    ### $f->load_id, "\n") if $verbose;
	    for my $mRNA ( @mrna ) {  
		my $cds;
		my ($id) = $mRNA->name || $mRNA->load_id;


		my $lastexon;
		my $i = 1;
		my @fs;
		if( $mRNA->has_tag('frameshifts') ) {
		    @fs = $mRNA->get_tag_values('frameshifts');
		}
		my $translation;
		if( $mRNA->has_tag('translation') ) {
		    ($translation) = $mRNA->get_tag_values('translation');
		}
		my $transl_table = 1;
		if( $mRNA->has_tag('transl_table') ) {
		    ($transl_table) = $mRNA->get_tag_values('transl_table');
		}
		my @CDS = sort { $a->strand * $a->start <=>
				     $b->strand * $b->start } $mRNA->get_SeqFeatures('cds');
		warn(" there are ", scalar @CDS, " cds features for $id\n");
		next if( @CDS == 0 );
		my $splitloc = Bio::Location::Split->new();
		my $singlexon = (@CDS > 1) ? 0 : 1;
		if( $singlexon && 
		    defined $CDS[0]->phase &&
		    $CDS[0]->phase == 1 ) {
		    if( $CDS[0]->strand > 0 ) {
			$CDS[0]->start( $CDS[0]->start + $CDS[0]->phase);
		    } else {
			$CDS[0]->end( $CDS[0]->end - $CDS[0]->phase);
		    }
		}
		for my $exon ( @CDS ) {
		    $splitloc->add_sub_Location($exon->location);
#		if( $exon->strand < 0 ) {
#		    $cds .= $db->segment($exon->seq_id,
#					 $exon->end => $exon->start)->seq;
#		} else {
		    if( $exon->length == 0 ) {
			warn("strand is ", $exon->strand, 
			     " for mRNA that is on ", $mRNA->strand, "\n");
		    } else {
			$cds .= $exon->dna;
		    }
#		}
		    unless( @fs ) {
			if( $lastexon ) {
			    my ($seqid,$s,$e) = ($exon->seq_id);
			    my $loc;
			    if( $exon->strand > 0 ) {			
				($s,$e) = ($lastexon->end+1, $exon->start-1);
				$loc = Bio::Location::Simple->new(-start => $s,
								  -end   => $e);
			    } else {
				($s,$e) = ($lastexon->start-1,$exon->end+1);
				$loc = Bio::Location::Simple->new(-start => $e,
								  -end   => $s,
								  -strand => $exon->strand
								  );    
			    }
			    my ($iseq) = $db->segment($seqid,$s,$e);
			    if( ! defined $iseq ) {
				warn("cannot get seq for $seqid $s..$e\n");
			    } else {
				$iseq = $iseq->seq->seq;
			    }

			    my $intron = Bio::PrimarySeq->new
				(-seq => $iseq,
				 -id => "$prefix:$id.i".$i++,
				 -desc => sprintf("%s %s:%s",
						  $gname,
						  $seqid,$loc->to_FTstring));
			    if( $intron->length < 4) {
				warn("intron ",$intron->length, " too short for: ", 
				     $intron->display_id,"\n");
			    } else {
				$iout->write_seq($intron);
			    }
			}
			$lastexon = $exon;
		    }
		}
		if( defined $cds ) {
		    my $cdslen = length($cds);

#		if( $singlexon && $cdslen % 3 != 0 ) 
#		{		    
#		    warn("dropping ",$cdslen % 3," bases\n $id");
#		    substr($cds,0,$cdslen % 3,'');
#		}
		    my $locstr = $splitloc->to_FTstring();
		    my $cdsseq = Bio::PrimarySeq->new(-seq => $cds,
						      -id  => "$prefix:$id",
						      -description => 
						      sprintf("gene=%s %s:%s",
							      $gname, 
							      $f->seq_id,
							      $locstr));	
		    my $pepseq;
		    if( $f->seq_id eq 'SPMICG' || $f->seq_id =~ /mito/i ) {    
			$pepseq = $cdsseq->translate(-codontable_id => 4);
		    } elsif( $f->seq_id =~ /^calb|clus|cpara|cdub|ctro/ ) {
			$pepseq = $cdsseq->translate(-codontable_id => 12);
		    } else { 
			$pepseq = $cdsseq->translate(-codontable_id => $transl_table);
		    }
		    my $peptide = $pepseq->seq;
		    substr($peptide,length($peptide)-1,1,'') 
			if( rindex($peptide,"*") == length($peptide)-1);
		    $pepseq->seq($peptide);
		    if( $peptide =~ /\*/ ) {
			warn( ">$id $gname ",$f->seq_id," $locstr\n","$peptide\n");
			warn( ">$id\_cds\n", $cds,"\n");
		    } elsif( defined ($translation) &&
			     $translation ne $peptide ) {
			warn("translation doesn't equal peptide:\n",
			     ">$id\_translation $locstr\n",$translation,"\n",
			     ">$id\_peptide $locstr\n",$peptide,"\n",
			     ">$id\_cds\n",$cds,"\n");
			$pepseq->seq($translation);
			$pepseq->description($cdsseq->description);
			$cout->write_seq($cdsseq);
			$pout->write_seq($pepseq);
		    } else {
			$pepseq->description($cdsseq->description);
			$cout->write_seq($cdsseq);
			$pout->write_seq($pepseq);
		    }
		} else {
		    warn("no defined cds for $id\n") if $verbose;
		}
	    }
	    last if $debug;
	}
    }
}
	 
