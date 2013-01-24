#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;

my $posfile = 'positive_polyA.fa';
my $negfile = 'negative_polyA.fa';
my $min_len_UTR = 100;
my $gff;
my $polyAfasta;
my $db;
my $polyA_window_width = 50;
my $upstream = 200;
my $downstream = 150;

GetOptions('p|pos:s' => \$posfile,
	   'n|neg:s' => \$negfile,
	   'g|gff:s' => \$gff,
	   'f|pasa|polyA:s' => \$polyAfasta,
	   'db:s'    => \$db,
	   'm|min:i' => \$min_len_UTR,
	   'w|window:i' => \$polyA_window_width,
	   'u|upstream:i' => \$upstream,
	   'd|downstream:i' => \$downstream,
    );

# arguments 
# genome
# GFF file with three_prime_utr features


if( ! defined $db ) {
    die("no DB provided\n");
}

# --[negative]--------AAAAA-----------[negative]
#   -200-150  (-50 -IN- +50) +150-200

my $dbh = Bio::DB::Fasta->new($db);
my $out_positive = Bio::SeqIO->new(-format => 'fasta',
				   -file   => ">$posfile");

my $out_negative = Bio::SeqIO->new(-format => 'fasta',
				   -file   => ">$negfile");
my @sites;
if( $gff ) {
    open(my $in => $gff) || die $!;
    my %count;
    while(<$in>) {
	next if /^\#/;
	next if /^\s+/;
	chomp;
	my ($seqid,$src,$type,$start,$end, $score, $strand,$frame,$group) = 
	    split(/\t/,$_);
	next unless $type eq 'three_prime_utr';
	next unless abs($end - $start) > $min_len_UTR;
	my %group = map { split (/=/,$_) } split(';',$group);
	my $id = $group{'Parent'};
	unless (defined $id ) { die "no Parent tag for feature $_\n"}
	if( $count{$id}++ ){
	    $id .= $count{$id};
	}
	push @sites, [$seqid,$start, $end, $strand, $id];
    }
} elsif( $polyAfasta ) {
    my $seqin = Bio::SeqIO->new(-format => 'fasta', -file => $polyAfasta);
    while( my $seq = $seqin->next_seq ) {
	my $id = $seq->id;
	my ($seqid,$start,$strand);
	if( $id =~ /^(\w+)\-(\d+)\_([\+\-])/ ) {
	    ($seqid,$start,$strand) = ($1,$2,$3);
	} else { 
	    warn("cannot match id $id\n");
	    next;
	}
	my $desc = $seq->desc;
	my @l = split(/\s+/,$desc);
	push @sites, [$seqid,$start,$start,$strand,pop @l];
    }
}

for my $site ( @sites ) {
    my ($seqid,$start,$end,$strand,$id) = @$site;
    my $seqlen = $dbh->length($seqid);

    my ($posseq,@negseqs);
    if( $strand eq '-' ) {
	my ($left,$right) = ( $start - $polyA_window_width,
			      $start + $polyA_window_width-1);
	$right = $seqlen if $right > $seqlen;
	$left = 1 if $left <= 0;	

	my $posseqstr = $dbh->seq($seqid, $left => $right);
	$posseq = Bio::Seq->new(-seq => $posseqstr,
				-id  => "$id.polyA",
				-desc => sprintf("%s:%d..%d strand=rev",
						 $seqid,$right, $left));
	$posseq = $posseq->revcom;
	$right = $start - $downstream;
	$left = $right - 2*$polyA_window_width+1;  
	my $negseqstr = $dbh->seq($seqid, $left => $right);
	my $negseq = Bio::Seq->new(-seq => $negseqstr,
				   -id  => "negDN_$id.polyA",
				   -desc => sprintf("%s:%d..%d strand=rev",
						    $seqid,$right, $left));
	push @negseqs, $negseq->revcom;

	$left = $start + $upstream;
	$right = $left + 2*$polyA_window_width-1;  
	$negseqstr = $dbh->seq($seqid, $left => $right);
	$negseq = Bio::Seq->new(-seq => $negseqstr,
				-id  => "negUP_$id.polyA",
				-desc => sprintf("%s:%d..%d strand=rev",
						 $seqid,$right, $left));
	push @negseqs, $negseq->revcom;
    } else {
	my ($left,$right) = ( $end - $polyA_window_width,
			      $end + $polyA_window_width-1);

	my $posseqstr = $dbh->seq($seqid, $left => $right);
	$posseq = Bio::Seq->new(-seq => $posseqstr,
				-id  => "$id.polyA",
				-desc => sprintf("%s:%d..%d strand=fwd",
						 $seqid,$left,$right));

	$left = $start - $upstream;
	$right = $left + 2*$polyA_window_width-1;  
	my $negseqstr = $dbh->seq($seqid, $left => $right);
	my $negseq = Bio::Seq->new(-seq => $negseqstr,
				   -id  => "negUP_$id.polyA",
				   -desc => sprintf("%s:%d..%d strand=fwd",
						    $seqid, $left,$right));
	push @negseqs, $negseq;

	$left = $start + $downstream;
	$right = $left + 2*$polyA_window_width-1;  
	$negseqstr = $dbh->seq($seqid, $left => $right);
	$negseq = Bio::Seq->new(-seq => $negseqstr,
				-id  => "negDN_$id.polyA",
				-desc => sprintf("%s:%d..%d strand=fwd",
						 $seqid, $left,$right));
	push @negseqs, $negseq;
    }
    $out_positive->write_seq($posseq);
    $out_negative->write_seq(@negseqs);
}
