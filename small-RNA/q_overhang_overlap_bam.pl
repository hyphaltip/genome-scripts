#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
use List::Util qw(sum);

use Env qw(HOME);
# provide SAM file(s) and look for evidence of siRNAs 
my (@bams,$genome);
my $debug = 0;
my $min_len = 21;
my $max_len = 28;
my ($start_off,$stop) = (308970,320950);
my $seqid = 'supercont10.6';
GetOptions(
    'v|verbose!'  => \$debug,    
    'b|bam:s'     => \@bams,
    'g|genome:s'  => \$genome,
    'start:i'     => \$start_off,
    'stop:i'      => \$stop,
    'seqid:s'     => \$seqid,
    );
for my $bam ( @bams ) {
    
    my (undef,undef,$bamshort) = File::Spec->splitpath($bam);
    my %offset;
    if( $bamshort =~ /(\S+)\.bam/ ) {
	my $stem = $1;
	$stem =~ s/\.clp_\d+\S*//;
	my $bam_db = Bio::DB::Sam->new(-bam => $bam,
				       -expand_flags => 1,
					-autoindex => 1);
				       #-fasta => $genome);
    
	my @alignments = $bam_db->get_features_by_location
	    (-seq_id => $seqid,
	     -start  => $start_off,
	     -end    => $stop );
	for my $aln ( @alignments ) {
    my %seen;
	    my $segment = $bam_db->segment($aln->seq_id,
					   $aln->start => $aln->end);
	    next if $aln->query->length > $max_len || $aln->query->length < $min_len;
	    my @qid = split(/:/,$aln->query->seq_id);
	    $qid[-1] =~ s/\#\d+$//;
	    shift @qid;
	    shift @qid;
	    my $qid = join(":",@qid);
	    next if $seen{$qid}++;
	    if( $debug ) {
		printf "q: %-20s %d..%d %d..%d %s%s (%s)\n",$qid,
		$aln->start, $aln->end,
		$aln->query->start, $aln->query->end,
		' 'x ($aln->start - $start_off),
		$aln->query->dna, $aln->query->strand;
	    }
	    for my $overlapread ( $segment->features ) {
	        next if $overlapread->query->strand == $aln->query->strand;
	        next if $overlapread->query->length > $max_len || $overlapread->query->length < $min_len;
		my @tid = split(/:/,$overlapread->query->seq_id);
		$tid[-1] =~ s/\#\d+$//;
		shift @tid;
		shift @tid;
		my $tid = join(":",@tid);
		next if $seen{$tid};
		my $space = ' 'x ($overlapread->start - $start_off);		
		
		if( $debug ) {
		    printf "t: %-20s %d..%d %d..%d %s%s (%s) [off: %d]\n", 
		    $tid,
		    $overlapread->start, $overlapread->end,
		    $overlapread->query->start, $overlapread->query->end,
		    $space, $overlapread->query->dna, $overlapread->query->strand,
		    $aln->start - $overlapread->start,
		    ;
		}
		$offset{$aln->start - $overlapread->start}++;
	    }
	    $start_off = $aln->start;
	    print "\n" if $debug

	}    
	open(my $fh => ">$stem.offsets") || die $!;
	for my $off ( sort { $a <=> $b } keys %offset ) {
	    print $fh join("\t", $off, $offset{$off}), "\n";
	}
    }
}
