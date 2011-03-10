#!usr/bin/perl -w
use strict;

=head1 NAME

search_gaps_for_TIRs - find based on TIR and TSD extraction

=head1 SYNOPSIS

 perl search_gaps_for_TIRs.pl gaps_file.fa [-max max_TIR_Length]

=head1 DESCRIPTION

 Provide a database of the query sequences 

=cut

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Perl;
use Bio::Tools::Run::WrapperBase;
use Getopt::Long;

my $Min_PID = 90; # require a minimum percent identity
my $TIR_max_search = 35;
my $TSD_max_search = 2;
my $TIR_min_length = 9;
my $TSD_f = 'TA';
my $TSD_r = reverse($TSD_f);
my $TSD_length = length($TSD_f);

my $tempfile_left = "/tmp/$$.left.fa";
my $tempfile_right = "/tmp/$$.right.fa";
my $wrapper = Bio::Tools::Run::WrapperBase->new;
my $bl2seq = $wrapper->executable('bl2seq'); # find bl2seq if it is in your path
my $noblast = 0;
my $db;
my $genome;
my $debug = 0;
GetOptions('d|db:s'       => \$db,
	   'g|genome:s'   => \$genome,
	   'n|noblast!'   => \$noblast,
	   'v|verbose!'   => \$debug,
);
$db ||= shift @ARGV;

my $dbh = Bio::DB::Fasta->new($genome);
my $seqio = Bio::SeqIO->new(-format => 'fasta', -file => $db);
my ($base) = $db;
$base =~ s/\.(fasta|fas|fa|seq|nt|dna)$//;
my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$base.TE_cands.fas");
my $ct;
while( my $gapseq = $seqio->next_seq ) {
    my ($lefthalf,$righthalf);
    if( $gapseq->length < 2*$TIR_max_search ) {
	$lefthalf = $gapseq->trunc(1,$TIR_max_search);
	$righthalf = $gapseq->trunc($gapseq->length - $TIR_max_search,$gapseq->length)->revcom;
    } else {
	$lefthalf = $gapseq->trunc(1,$TIR_max_search);
	$righthalf = $gapseq->trunc($gapseq->length - $TIR_max_search,$gapseq->length)->revcom;
    }
    my ($chrom,$loc) = split(/:/,$gapseq->description);
    my ($start,$end) = split('\.\.',$loc);
    
    $lefthalf->id('left');
    $righthalf->id('right');
    next if ($lefthalf->seq =~ tr/N/N/ > 3);

#    my ($f1,$tf1) = $wrapper->io->tempfile(CLEANUP=>1);
#    my ($f2,$tf2) = $wrapper->io->tempfile(CLEANUP=>1);
    
    my ($longest_com_sub,$length_TIR,$l_pos,$r_pos);
    Bio::SeqIO->new(-file => ">$tempfile_left",-format=>'fasta')->write_seq($lefthalf);
    Bio::SeqIO->new(-file => ">$tempfile_right",-format=>'fasta')->write_seq($righthalf);    
#    close($f1);
#    close($f2);
    open(my $run => "bl2seq -i $tempfile_left -j $tempfile_right -p blastn -D 1 -F F|");
    my @best;
    # save the longest alignment;
    while (<$run>) {
	next if /^\#/;
	my @row = split;
	next if $row[2] < $Min_PID; # require a minimum percent ID
	if ( ! defined $best[0] || $best[3] < $row[3] ) {
	    @best = @row;
	}
    }
    
    if ( @best && !( $best[6] == 1 && $best[7] >= $lefthalf->length ) ) {
	$length_TIR = $best[3];
	$l_pos = $best[6];
	$r_pos = $best[8];
	$longest_com_sub = $lefthalf->subseq($best[6],$best[7]);
    } else { next }
    my ($TSD_search_left) = $dbh->seq($chrom,$start - $TSD_max_search+1,$start);
    my ($TSD_search_right) = $dbh->seq($chrom,$end, $end +$TSD_max_search-1);
    
    next if $longest_com_sub && lc($longest_com_sub) eq $longest_com_sub;
    next unless $length_TIR && $l_pos && $r_pos;
    
    if ( $length_TIR > $TIR_min_length ) { #&&
	my $TIR_right = ( $best[8] > $best[9] ) ? $righthalf->trunc(sort { $a <=> $b } $best[8],$best[9])->revcom->seq : $righthalf->subseq($best[8],$best[9]);

      my $TIR_left = $lefthalf->subseq($best[6],$best[7]);
      

      my $TR_type = $best[6] < $best[7] && $best[8] < $best[9] ? 'TIR' : 'TDR';
      print join("\t", $gapseq->id, $chrom, "$start..$end", $gapseq->length, 
		 $length_TIR, $TR_type, $TIR_left,$TIR_right,$TSD_search_left, $TSD_search_right, 
		 $TSD_search_left eq $TSD_search_right ? 'TSD_MATCH' : 'TSD_NOMATCH'),"\n";
	my $composite_seq = uc($TSD_search_left) . lc($TIR_left) . 
	    $dbh->seq($chrom,$start+$length_TIR+1,$end- $length_TIR-1).
	    lc (revcom_as_string($TIR_right)) .
	    uc($TSD_search_left);
	my $desc= sprintf("%s:%d..%d",$chrom,
			  $start - length($TSD_search_left),
			  $end + length($TSD_search_right));
	$out->write_seq(Bio::Seq->new(-seq  => $composite_seq,
				      -id   => "$base\_TE_CAND_".$ct++,
				      -desc => $desc));
  }
    last if $debug;
}

warn "$base -- $ct Indels have TIRs\n";

