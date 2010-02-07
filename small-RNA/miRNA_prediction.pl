#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::Sam;
use Bio::DB::Fasta;
use Bio::SeqIO;
use List::Util qw(sum);
use IPC::Open2;
use File::Spec;
use RNA;

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
my $min_coverage = 20;
my ($window_table,$bam,$fa);
my ($upstream,$downstream) = (160,20);
my $target_dir = 'targets';
my $min_mfe = -18;
my $min_stemlen = 18;
my $debug = 0;
GetOptions(
	   'v|verbose!'    => \$debug,
	   'w|t|table:s'   => \$window_table,
	   'b|bam:s'       => \$bam,
	   'fa|g|genome:s' => \$fa,
	   'm|min:i'       => \$min_coverage,
	   'mfe:s'         => \$min_mfe,
	   'stem:i'        => \$min_stemlen,
	   );

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
while(<$fh>) {
    my ($seqid,$start,$end,@window_counts) = split;
    next unless ( sum(@window_counts) > $min_coverage );
    $count++;
    
    for my $aln ( $db->get_features_by_location
		  (-seq_id => $seqid,
		   -start  => $start,
		   -end    => $end)) 
    { 
	my $aln_count = 1;
	if( $aln->name =~ /^(\d+)\-(\d+)$/) {
	    $aln_count = $2;
	}    
	my $astart = $aln->start;
	my $aend   = $aln->end;
	my $pred_segment_fwd = $dbfa->seq($seqid,
					    $astart - $upstream,
					    $aend + $downstream);
	
	my $pred_segment_rev = $dbfa->seq($seqid,
					  $astart - $downstream,
					  $aend + $upstream);
	for my $coords ( [$astart-$upstream, $aend], 
			 [$aend+$upstream,   $astart-$downstream] ) {
	    my $seq = $dbfa->seq($seqid,@$coords);

	    my ($struct, $mfe) = RNA::fold($seq); 
	    my $name_window = "$seqid:".join("..",@$coords);
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
		    my ($e,$str,$prob,$shape) = split(/\s+/,$_,4);
		    my @shapes;
		    my $stemlen;
		    if( $shape eq '[]' ) {
			if ($str =~ /^\.*([\(\.]+\()(\.+)(\)[\.\)]+\))\.*/) {
			    my ($steml,$loop,$stemr) = ($1,$2,$3);
			    warn "$steml   $loop   $stemr for $str\n" if $debug;
			    $stemlen =  max(length($steml),length($stemr));
			} else {
			    warn("no stem parsed in $str\n");
			}
			
			@shapes = (#sprintf("ShapesEnergy=%.2f",$e),
				   # sprintf("ShapesStructure=%s",$str),
				   sprintf("Shape=%s",$shape),
				   sprintf("ShapesProb=%.5f",$prob),
				   sprintf("StemLen=%d",$stemlen),
				   );
			last;
		    }
		    print $ofh join("\t", 
				    $seqid,'miRNAPred','ncRNA',
				    (sort { $a <=> $b} @$coords),$prob,
				    '.', '.',
				    join(";",@shapes)),"\n";
		    next unless $stemlen > $min_stemlen;
		    close($rdr);
		    close($wtr);
		    my $systematic_name = sprintf($pref,$mirnacounter++);
		    my $seqobj = Bio::Seq->new(-id   => $systematic_name,
					       -desc => sprintf($name_window),
					       -seq => $seq);
		    $outseq->write_seq($seqobj);	
		    waitpid $pid, 0;
		    exit;
		}
	    }		
	}
	last if $debug;
    }
    last if $debug;
    warn($count, " windows\n");
}
print "$count windows\n";
