#!/usr/bin/perl -w
use strict;

use Env;
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
my $longest;
my $user = $USER;
my $pass;
my $dointrons = 0;
my ($verbose,$force) = (0,0);
GetOptions(
	   'f|force!'   => \$force,
	   'd|dir:s'    => \$dir,
	   'e|ext:s'    => \$ext,
	   'c|cds:s'    => \$cdsdir,
	   'p|pep:s'    => \$pepdir,
	   'i|intron:s' => \$introndir,
	   'g|gene:s'   => \$genedir,
	   'v|verbose'  => \$verbose,
	   'u|user:s'        => \$user,
	   'pass|password:s' => \$pass,
	   
	   'l|longest!'  => \$longest,
	   'di|dointrons!'  => \$dointrons,
	   );
$dir = shift @ARGV unless $dir;

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

    my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
					     -dsn     => "dbi:mysql:$genome",
					     -user    => $user,
					     -pass    => $pass);
    my $pout = Bio::SeqIO->new(-format => 'fasta',
			       -file   => ">$pepdir/$genome.pep.fa");
    
    my $cout = Bio::SeqIO->new(-format => 'fasta',
			       -file   => ">$cdsdir/$genome.cds.fa");
    
    my $iout = Bio::SeqIO->new(-format => 'fasta',
			       -file   => ">$introndir/$genome.intron.fa");

    my $gout = Bio::SeqIO->new(-format => 'fasta',
			       -file   => ">$genedir/$genome.gene.fa");
    
    my @features = $db->get_features_by_type('gene');
    warn("there are ", scalar @features, " gene features\n") if $verbose;
    for my $f ( @features ) {
	my ($gname) = $f->load_id;
	my @mrna = $f->get_SeqFeatures('mRNA');
	### warn("there are ", scalar @mrna, " features for gene ", 
	### $f->load_id, "\n") if $verbose;
	my %longest_transcript;
	
	for my $mRNA ( @mrna ) {  
	    my $cds;
	    my ($id) = $mRNA->load_id;
	    
	    
	    my $lastexon;
	    my $i = 1;
	    my @fs;
	    my $translation;
	    if( $mRNA->has_tag('translation') ) {
		($translation) = $mRNA->get_tag_values('translation');
	    }
	    if( $mRNA->has_tag('frameshifts') ) {
		@fs = $mRNA->get_tag_values('frameshifts');
	    }
	    
	    my @CDS = $mRNA->get_SeqFeatures('cds');
	    my @introns;
	    ### warn(" there are ", scalar @CDS, " cds features for $id\n");
	    my $splitloc = Bio::Location::Split->new();
	    for my $exon ( sort { ($a->start*$a->strand) <=> 
				      ($b->end*$b->strand) } 
			   @CDS ) {		
		if( $exon->length == 0 ) {
		    warn("strand is ", 
			 $exon->strand, " for mRNA that is on ", 
			 $mRNA->strand, "\n");
		    
		} else {
		    $cds .= $exon->dna;
		}
		$splitloc->add_sub_Location($exon->location);
		if( ! @fs && $dointrons ) {
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
			my $iseq = $db->segment($seqid,$s,$e);
			if( ! $iseq ) {
			    warn("cannot get seq for $seqid\n");
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
			    push @introns, $intron;
			    # $iout->write_seq($intron);
			}
		    }
		    $lastexon = $exon;
		}
	    }
	    if( defined $cds ) {
		my $locstr = $splitloc->to_FTstring();
		my $cdsseq = Bio::PrimarySeq->new(-seq => $cds,
						  -id  => "$prefix:$id",
						  -description => 
						  sprintf("%s %s:%s",
							  $gname, 
							  $f->seq_id,
							  $locstr));	
		my $pepseq;
		if( $f->seq_id eq 'SPMICG' || $f->seq_id =~ /mito/i || $f->seq_id =~ /^MT/i) {    
		    $pepseq = $cdsseq->translate(-codontable_id => 4);
		} elsif( $f->seq_id =~ /^calb|clus|cpara|cdub|ctro/ ) {
		    $pepseq = $cdsseq->translate(-codontable_id => 12);
		} else { 
		    $pepseq = $cdsseq->translate;
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
		    
		} else {
		    my $total_length = $pepseq->length;
		    $pepseq->description($cdsseq->description);
		    if( ! defined $longest_transcript{'length'} ||
			$longest_transcript{'length'} < $total_length ) {
			$longest_transcript{'length'} = $total_length;

			$longest_transcript{'cds'} = $cdsseq;
			$longest_transcript{'pep'} = $pepseq;
			
			$longest_transcript{'introns'} = [@introns];
		    }
#		    $cout->write_seq($cdsseq);
#		    $pout->write_seq($pepseq);
		}
	    } else {
		warn("no defined cds for $id\n") if $verbose;
	    }
	}
	if( defined $longest_transcript{'length'} ) {
	    $cout->write_seq($longest_transcript{'cds'});
	    $pout->write_seq($longest_transcript{'pep'});
	    $iout->write_seq(@{$longest_transcript{'introns'} || []});
	}
    }
}
	 
