#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;

my $odir = 'orthologs';
my $db;
my $expected_spcount = 13;
GetOptions('o|out:s' => \$odir,
	   'db:s'    => \$db,
	   'c|count:i'=> \$expected_spcount,
	   );

warn("odir is $odir\n");
unless( defined $db ) {
    die("need to provide a db with --db\n");
}

mkdir($odir) unless -d $odir;
my $dbh = Bio::DB::Fasta->new($db);

while(<>) {
    my ($group,@orthologs) = split;
    $group =~ s/://;
    my %count;
    for my $gene ( @orthologs ) {
	my ($sp)= split(/\|/,$gene);
	$count{$sp}++;
    }
    if( keys %count == $expected_spcount ) { 
	my $only_1 = 1;
	for my $ct ( values %count ) {
	    if( $ct != 1 ) { # see if every sp only has 1 gene
		$only_1 = 0;
		last;
	    }
	}
	if( $only_1 ) {
	    my $out = Bio::SeqIO->new(-format => 'fasta',
				      -file   => ">$odir/$group.fa");
	    for my $gene ( @orthologs ) {
		my ($sp,$id) = split(/\|/,$gene);
		my $seq = Bio::Seq->new(-seq => $dbh->seq($gene),
					-id  => $sp,
					-desc=> $id);
		$out->write_seq($seq);
	    }
	}
    }
}
