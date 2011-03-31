#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::PrimarySeq;
my $db;
my $col_count = 6;
my $outdir = 'merge_orthologs';
GetOptions(
	   'd|dir|db:s' => \$db,
	   'o|out:s'    => \$outdir,
	   );

unless( $db && (-d $db || -f $db) ) {
    die("cannot open a db or db ($db) directory for CDS sequence");
}
mkdir($outdir);
my $dbh = Bio::DB::Fasta->new($db);
my @header = split(/\s+/,<>);

while(<>) {
    my @row = split;
    my @group = splice(@row,0,$col_count);
    my ($gene_id,$ref) = @group;
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$outdir/$gene_id.cds.fa");
    my @header_r = @header;
    my @h = splice(@header_r,0,$col_count);
    my $cds = $dbh->seq($ref);
    my $seq = Bio::PrimarySeq->new(-id => 'ncra10',
				   -description => 
				   sprintf("gene=%s mRNA=%s %s:%d..%d %s",
					   @group),
				   -seq => $cds);
    $out->write_seq($seq);
    while( @row ) {
	@group = splice(@row,0,$col_count);
	@h = splice(@header_r,0,$col_count);
	my ($pref) = split(/_/,$h[0]); # species prefix is first bit before _
	my $mRNA = $group[1]; # 2nd col in a group is always the mRNA
	# which is what we named all the transcript-protein IDs
	$cds = $dbh->seq("$pref:$mRNA");
	$seq = Bio::PrimarySeq->new(-id => $pref,
				    -description => 
				    sprintf("gene=%s mRNA=%s %s:%d..%d %s",
					    @group),
				    -seq => $cds);
	
	if( ! defined $cds ) {
	    warn("cannot find cds for $pref:$mRNA ($pref)");
	    next;
	}
	$out->write_seq($seq);
    }
}



