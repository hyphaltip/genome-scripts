#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::Sam;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Collection;
use List::Util qw(sum max);
use IPC::Open2;
use File::Spec;
use RNA;

# 
# At each genomic locus of a short sequence to be tested, two
# sequences covering the read were extracted for secondary structure
# analysis, one sequence extending 160 nt upstream and 20 nt
# downstream from the read, and the other covering 20 nt upstream and
# 160 nt downstream of the read. The secondary structures for these
# two sequences were predicted by the RNA-fold program [61]. Those
# sequences that met the following two criteria were considered as
# candidate miRNA precursors. First, a secondary structure must have a
# hairpin with at least 18 paired nucleotides in its stem
# region. Second, the hairpin must have free energy less than or equal
# to -18 kCal/mol. and one central loop. In the biogenesis of miRNAs,
# other portions of pre-miRNAs will degrade soon after DCL1 separates
# miRNA:miRNA* duplex from the pre-miRNAs. Thus, the reads detected in
# the 454 libraries are supposed to be mature miRNAs or miRNA*. We
# thus checked the reads that are supposed to be mature miRNA or
# miRNA* to satisfy the requirements of the MIRCHECK program [31].

my $pref = 'sRNAPred%04d';
my $Min_coverage = 50;
my ($window_table,$bam,$fa);
my ($upstream,$downstream) = (160,20);
my $target_dir = 'targets';
my $min_mfe = -40;
my $min_stemlen = 18;
my $debug = 0;
GetOptions(
	   'v|verbose!'    => \$debug,
	   'w|t|table:s'   => \$window_table,
	   'b|bam:s'       => \$bam,
	   'fa|d|db|g|genome:s' => \$fa,
	   'min:i'         => \$Min_coverage,
	   'mfe:i'         => \$min_mfe,
	   'stem:i'        => \$min_stemlen,
	   );
#if( $min_mfe > 0 ) { 
#    $min_mfe = "-".$mfe;
#}
my $max_tile = $upstream + $downstream + 30;

$bam ||= shift @ARGV;
$fa  ||= shift @ARGV;

if( ! defined $window_table ) {
    warn("cannot proceed without a window_table -w/-t\n");
    exit;
}
my $dbfa = Bio::DB::Fasta->new($fa);
my $db = Bio::DB::Sam->new(-bam => $bam,
			   -fa  => $fa);

my ($bamstem) = ( $bam =~ /(\S+)\.bam$/);
$bamstem =~ s/\.clp_collap//;
my $outseq = Bio::SeqIO->new(-format => 'fasta',
			     -file   => ">$bamstem.miRNApred.fas");
open(my $ofh => ">$bamstem.miRNApred.gff3") || die $!;
open(my $fh => $window_table ) || die $!;
my $header = <$fh>;
my @header = split(/\t/,$header);

my $mirnacounter = 0;
my $count = 0;
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -file   => ">$bamstem.folding_seqs.fas");
while(<$fh>) {
    my ($seqid,$start,$end,@window_counts) = split;    
    next unless ( &min_coverage(@window_counts) );    
    warn("here with $seqid $start $end @window_counts\n") if $debug;
    $count++;
    
    my @tile = &tile_features($db->get_features_by_location
			      (-seq_id => $seqid,
			       -start  => $start,
			       -end    => $end));
    for my $tile_window ( @tile ) {
	my ($astart,$aend,$aln_count) = @$tile_window;
	my ($pred_segment_fwd, $pred_segment_rev);
	my @coords;
	if( abs($astart - $aend) > $max_tile ) {
	    push @coords, ([$astart, $aend,$aln_count], [$aend, $astart,$aln_count]);
	} else {
	    push @coords, ( [$astart - $upstream,
			     $aend + $downstream,
			     $aln_count], 
			    [$astart - $downstream,
			    $aend + $upstream,
			     $aln_count]);
	}
	
	for my $coords ( @coords ) {
	    my ($s,$e,$ct) = @$coords;
	    my $seq = $dbfa->seq($seqid,$s => $e);
	    my ($struct, $mfe) = RNA::fold($seq); 
	    my $name_window = "$seqid\_".join("-",$s,$e);
	    $out->write_seq(Bio::Seq->new(-id => $name_window,
					  -desc=> "DEPTH=$ct MFE=$mfe",
					  -seq => $seq));
	    
	    if( $mfe < $min_mfe )  {
#		if( $struct =~ /([(.]+)(\.{2,})([.)]+)/ ) {		
#		    my ($left,$loop,$right)= ($1,$2,$3);
#		    my $leftct = ($left =~ tr/\(/\(/);
#		    if( $leftct > $min_stemlen ) {			
#			warn("$leftct $struct, $mfe\n");
#		    }
#			    
		my ($rdr,$wtr);
		my $pid = open2($rdr,$wtr,"RNAshapes", "-q","-m","[]", "-F","0.1");
		Bio::SeqIO->new(-format => 'fasta',
				-fh     => $wtr)->write_seq
				    (Bio::Seq->new(-seq => $seq,
						   -id  => $name_window));
		while(<$rdr>) {
		    next if( /^>/ || /^\s+(\S+)\s+/ );
		    next if /shape\s+\S+\s+not found/;
		    chomp;
		    my ($energy,$str,$prob,$shape) = split(/\s+/,$_,4);
		    my @shapes;
		    my $stemlen;
		    if( $shape eq '[]' ) {
			if ($str =~ /^\.*([\(\.]+\()(\.+)(\)[\.\)]+\))\.*/) {
			    my ($steml,$loop,$stemr) = ($1,$2,$3);
			    #warn "$steml   $loop   $stemr for $str\n" if $debug;
			    $stemlen =  max(length($steml),length($stemr));
			} else {
			    warn("no stem parsed in $str\n");
			}
			
			@shapes = (#sprintf("ShapesEnergy=%.2f",$e),
				   # sprintf("ShapesStructure=%s",$str),
				   sprintf("Shape=%s",$shape),
				   sprintf("Depth=%s",$ct),
				   sprintf("ShapesProb=%.5f",$prob),
				   sprintf("Energy=%s",$energy),
				   sprintf("StemLen=%d",$stemlen),
				   );
			my $systematic_name = sprintf($pref,$mirnacounter++);
#		    next unless $stemlen > $min_stemlen;
			print $ofh join("\t", 
					$seqid,'miRNAPred','ncRNA',
					(sort { $a <=> $b} ($s,$e)),$prob,
					'.', '.',
					join(";",sprintf("ID=%s;Name=%s",$systematic_name,$systematic_name),@shapes)),"\n";
			close($rdr);
			close($wtr);
			my $seqobj = Bio::Seq->new(-id   => $systematic_name,
						   -desc => sprintf("LOC=%s:%d..%d DEPTH=%d MFE=%s SHAPE_MFE=%s SHAPE_PROB=%s", 
								    $seqid,$s,$e,$ct,$mfe,$energy,$prob),
						   -seq => $seq);
			$outseq->write_seq($seqobj);	
			waitpid $pid, 0;
		    }
		}
	    }
	}
	last if $debug && $count > 8;
    }
    last if $debug && $count > 8;
   
}
print "$count windows\n";

sub min_coverage {
    for my $r ( @_) { 
	if( $r > $Min_coverage ) {
	    return 1;
	}
    }
    0;
}

sub tile_features {    
    my @features = @_;
    my %clusterlookup;
    my $clusterct;
    my @clusters;
    for my $aln1 ( @features ) {
	my $name1 = $aln1->name;
	my $bin = $clusterlookup{$name1};
	unless( defined $bin ) {
	    $bin = $clusterlookup{$name1} = $clusterct++;
	    push @{$clusters[$bin]}, $aln1;
	}
	for my $aln2 ( @features ) {
	    my $name2 = $aln2->name;
	    next if $name1 eq $name2;
	    next unless &aln_overlaps($aln1,$aln2);	    
	    my $bin2 = $clusterlookup{$name2};
	    if( defined $bin2 ) {
		if( $bin != $bin2 ) {
		    # join two bins
		    for my $rename_aln ( @{$clusters[$bin2]} ) {
			$clusterlookup{$rename_aln->name} = $bin;
		    }
		    push @{$clusters[$bin]}, @{$clusters[$bin2]};
		    $clusters[$bin2] = [];
		} 
		next;		
	    } else {
		push @{$clusters[$bin]}, $aln2;
	    }
	    $clusterlookup{$name2} = $bin2 = $bin;
	}
    }
    my @final_clusters;
    for my $cl (  grep { defined && scalar @$_ > 0} @clusters) {
	my ($start,$end,$incount);
	for my $aln ( @$cl ) {
	    $start = $aln->start unless defined $start && $start < $aln->start;
	    $end = $aln->end unless defined $end && $end < $aln->end;
	    if( $aln->name =~ /^(\d+)\-(\d+)$/) {
		$incount += $2;	    
	    } else { $incount++ }
	}
	push @final_clusters, [$start,$end,$incount];
    }
    return sort { $b->[2] <=> $a->[2] } @final_clusters;
}

sub aln_overlaps {
    my ($aln_1, $aln_2) = @_;
    
    return not ( $_[0]->start > $_[1]->end ||
		 $_[0]->end   < $_[1]->start);
}
