#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

use constant CDSEXON => 'cds';
use constant EXON => 'exon';
use constant MRNA => 'mRNA';
use constant GENE => 'gene';
my $verbose = 0;

my $write_sequence;
GetOptions("w|s|write|seq!" => \$write_sequence,
	   'v|verbose!'     => \$verbose,
	   );

my $gffdir = "raw/jgi";
my $ntdir = "raw/nt";

my $odir = "genomes";
my $src_string = "JGI";
mkdir($odir) unless -d $odir;

opendir(DIR, $gffdir) || die $!;

FILE: for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)_transcripts\.g[tf]f$/;
    my $stem = $1;
    warn("file is $file\n") if $verbose;

    my ($genus,$species,$strain,$version) = split(/_/,$stem);
	
    if( ! defined $species ) { 
	warn("cannot parse genus-species from $stem\n");
	next;
    }
    $version = $strain if ! defined $version;
    
    if( ! defined $version ) {
	warn("cannot find version in $stem\n");
	next;
    }

    mkdir("$odir/$stem") unless -d "$odir/$stem";
    if(  -f "$odir/$stem/$stem\_transcripts.gff3" &&
	 ! -z "$odir/$stem/$stem\_transcripts.gff3") {
	warn(" skipping $stem - already processed $stem\_transcripts.gff3\n");
	next;
    }
    
    my ($prefix) = lc substr($genus, 0, 1).substr($species,0,3);
    my $fafile;
    for my $bck ( qw(_supercontigs _contigs), '') {
	my $test = "$ntdir/$stem$bck.fasta";
	if( -f $test ) {
	    $fafile = $test;
	    last;
	}
    }
    if( ! defined $fafile ) {
	warn("cannot find a fasta file for $stem, skipping\n");
	next;
    }
    warn("fafile is $fafile\n") if $verbose;
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
				-file   => $fafile);
    my @seqs;
    my %seqids;
    
    while( my $seq = $seqio->next_seq ) {
	my $seqid = $seq->display_id;
	my $id;
	if( $seqid =~ /scaffold_(\d+)/ ) {
	    $id = sprintf("%s_%s",$prefix,$1);
        } elsif( $seqid =~ /^chr/ ) {
           $id = sprintf("%s_%s", $prefix,$seqid);
	} else {
	    warn("cannot fix ID for seqid: $seqid ($file)\n");
	    next FILE;
	}
	$seqids{$id} = $seq->length;
	$seq->display_id($id);
	$seq->description("");
	push @seqs, $seq;
    }

    my ($out,$in);
    open($in, "$gffdir/$file") || die $!;
    my (%genes);
    my %lookup;
    while(<$in>) { 
	my ($seqid,$src,$type,$start,$end,$score,
	    $strand,$frame,$lastcol) = split(/\t/,$_,9);
	next unless $type eq 'CDS' || $type eq 'exon';
    
	my $id;
	if( $seqid =~ /scaffold_(\d+)/ ) {
	    $id = sprintf("%s_%s",$prefix,$1);
        } elsif( $seqid =~ /^chr/ ) {
           $id  = sprintf("%s_%s", $prefix,$seqid);
	} else {
	    warn("cannot match contig out of $seqid\n");
	    close($in);
	    next FILE;
	}
	$seqid = $id;
	if( ! defined $seqids{$seqid} ) {
	    warn("no sequence for $seqid\n");
	    exit;
	} elsif ( $seqids{$seqid} < $end ) {
	    warn("seq length of $seqid is $seqids{$seqid} is less than annotation end ($end)\n");
	    exit;
	}
	
	my ($geneid,$transcriptid);
	if( $lastcol =~ /(?:name|gene_id)\s+\"([^\"]+)\";/ ) {
	    $geneid = $1;
	}
	
	if( $lastcol =~ /featureId\s+(\d+)/ ) {
	    $transcriptid = "$prefix\_t_$1";
	}
	if( $lastcol =~ /transcriptId\s+(\d+)/ ) {
	    $transcriptid = "$prefix\_t_$1";
	} elsif( $lastcol =~ /transcript_id\s+\"([^\"]+)\"/ ) { #"
	    $transcriptid = "$prefix\_t_$1";
        }
	
	if( defined $transcriptid && ! defined $lookup{$geneid}) {
	    $lookup{$geneid} = $transcriptid;
	}
	if( ! defined $transcriptid ) {
	    $transcriptid = $lookup{$geneid};
	}
	if( ! defined $geneid || ! defined $transcriptid ) {
	    warn("cannot find geneid or transcriptid in $lastcol\n");
	    next FILE;
	}
	my $typetag = $type eq 'CDS' ? CDSEXON : EXON;

	push @{$genes{$geneid}->{$transcriptid}}, 
	[$seqid,$typetag, $start,$end,$score,$strand,$frame];
    }

    mkdir("$odir/$stem") unless -d "$odir/$stem";
    open($out, ">$odir/$stem/$stem\_transcripts.gff3") || die $!;
    print $out "##gff-version 3\n";
    print $out "#date-prepared ".localtime(time)."\n";
    if( $write_sequence) {
	for my $seq ( @seqs ) {
	    print $out join("\t", $seq->display_id,
			    'chromosome',
			    'scaffold',
			    1, $seq->length,
			    '.', '+','.', sprintf("ID=%s;Name=%s",
						  $seq->display_id,$seq->display_id)),"\n";
	}
    }
    my %num= ('exon' => 1, 'cds' => 1);
    for my $gene ( keys %genes ) {
	my ($gmin, $gmax,$gstrand,$seqid_all);
	my @lines;
	for my $transcript ( keys %{$genes{$gene}} ) {
	    my ($tmin,$tmax,$tstrand);
	    for my $CDS ( @{$genes{$gene}->{$transcript}} ) {
		my ($seqid,$type,$start,$end,$score,$strand,$frame) = @$CDS;
		($start,$end) = sort { $a <=> $b } ($start,$end);
		$gmin = $start if ! defined $gmin || $gmin > $start;
		$gmax = $end if ! defined $gmax || $gmax < $end;
		$gstrand = $strand;

		$tmin = $start if ! defined $tmin || $tmin > $start;
		$tmax = $end if ! defined $tmax || $tmax < $end;
		$tstrand = $strand;
		
		push @lines, join("\t", $seqid, $src_string, $type,
				  $start, $end, $score, $strand, $frame,
				  join(';',
				       sprintf("ID=%s_%s%05d",$prefix,lc($type),
					       $num{$type}++),
				       sprintf('Parent=%s',$transcript)));
		$seqid_all = $seqid;
	    }
	    unshift @lines, join("\t", $seqid_all, $src_string, MRNA,
				 $tmin, $tmax, '.', $tstrand, '.',
				 join(';',
				      sprintf("ID=%s;Name=%s",$transcript,$transcript),
				      sprintf('Parent=%s',$gene)));
	}
	print $out join("\t", $seqid_all, $src_string, GENE,
			$gmin, $gmax, '.', $gstrand, '.',
			join(';',
			     sprintf("ID=%s",$gene),
			     sprintf("Name=%s.%d",$gene,$version))
			), "\n", join("\n", @lines), "\n";
    }
    if( $write_sequence ) {
	my $outio = Bio::SeqIO->new(-format => 'fasta',
				    -fh       => $out);
	$outio->write_seq(@seqs);
    }
}

