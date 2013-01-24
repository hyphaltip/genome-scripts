#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;

my %est;

my ($trimpoly,$gff_aln,$db,$genome) = @ARGV;
my ($queryfolder,$hitfolder,$alnfolder) = ('query','target','alns');

mkdir($queryfolder);
mkdir($hitfolder);
mkdir($alnfolder);
my $genome_dbh = Bio::DB::Fasta->new($genome);
my $dbh = Bio::DB::Fasta->new($db);

open(my $in => $trimpoly) || die $!;
my $out = Bio::SeqIO->new(-format => 'fasta',
			  -file   => ">$trimpoly.seqs");
while(<$in>) {
    my ($seqid,$perc_N,$end5,$end3,$initial_length) = split;
    my ($len) = abs($end3 - $end5) + 1;
    if( $initial_length != $len ) {
	$est{$seqid} = [$end5,$end3,$initial_length];    
#	print ("$seqid $end5 $end3 $initial_length\n");
	my @polyA;
	if( $end5 > 1 ) {
	    push @polyA, "5prime";
	}
	if( $end3 < $initial_length ) {
	    push @polyA, "3prime";
	}
	my $seqstr = $dbh->seq($seqid,$end5 => $end3);
	my $s = Bio::Seq->new(-seq => $seqstr,
				      -id  => $seqid,
			      -desc => join(",",@polyA));
	$out->write_seq($s);	
	Bio::SeqIO->new(-format => 'fasta',
			-file   => ">$queryfolder/$seqid")->write_seq($s);	
	Bio::SeqIO->new(-format => 'fasta',
			-file   => ">$queryfolder/$seqid.untrim")->write_seq(
			    $dbh->get_Seq_by_acc($seqid));	
    }
}

open($in => $gff_aln ) || die $!;
open(my $jobfh => ">run_exonerate.sh");
print $jobfh "module load exonerate\n";
my %seen;
while(<$in>) {
    next if /^\s+$/ || /^\#/;
    chomp;
    my ($seqid,$source,$type,$start,$end,$score,$strand,$frame,$group) = split(/\t/,$_);
    my %group = map { split(/=/,$_) } split(/;/,$group);
    my $target_str = $group{'Target'} || die "no target in feature\n";
    my ($target,$tstart,$tend,$tstrand) = split(/\s+/,$target_str);
    if( $est{$target} ) {
	printf "%s -> %d..%d (%d bp) to aln parts %d..%d (%s) on genome %s:%d..%d (%s)\n",
	$target,@{$est{$target}},
	$tstart,$tend,$tstrand,
	$seqid,$start,$end,$strand;
	if( ! -f "$hitfolder/$seqid" ) { 
	    my $outquery = Bio::SeqIO->new(-format => 'fasta',
					   -file   => ">$hitfolder/$seqid");
	    $outquery->write_seq($genome_dbh->get_Seq_by_acc($seqid));
	}	
	if( ! $seen{$target}++ ) {
	    #only need to request the alignment once
	    print $jobfh "exonerate -m est2genome --refine region --bestn 1 -q $queryfolder/$target -t $hitfolder/$seqid --showtargetgff > $alnfolder/$target.exonerate\n";
	}
    }    
    
}
