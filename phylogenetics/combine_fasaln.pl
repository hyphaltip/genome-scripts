#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Getopt::Long;
my $iformat = 'fasta';
my $oformat = 'nexus';
my $outfile = 'allseq.nex';
my $ext = 'fasaln.trim';
my $dir;
my @expected;
my $expected_file;
GetOptions('d|dir:s'   => \$dir,
	   'ext:s'     => \$ext,
	   'if:s'       => \$iformat,
	   'of:s'       => \$oformat,
	   'expected:s' => \$expected_file,
	   'o|out:s'   => \$outfile,
	   );

die("need a dir") unless $dir && -d $dir;

opendir(DIR, $dir) || die"$dir: $!";

if( $expected_file && open(my $fh => $expected_file) ) {
    while(<$fh>) {
	chomp;
	push @expected, $_;
    }
}
my (%matrix);

for my $file (sort readdir(DIR) ) {
    next if $file eq $outfile;
    next unless ($file =~ /(\S+)\.\Q$ext\E$/);
    my $in = Bio::AlignIO->new(-format => $iformat,
			       -file   => "$dir/$file");
    if( my $aln = $in->next_aln ) {
	my %seen;
	for my $seq ( $aln->each_seq ) {
	    my $id = $seq->id;
	    if( $id =~ /(\S+)\|/) { 
		$id = $1;
	    }
	    $matrix{$id} .= $seq->seq;
	    $seen{$id}++;
	}
	for my $exp ( @expected ) {
	    if( ! $seen{$exp} ) {
		$matrix{$exp} .= '-' x $aln->length;
	    }
	}
    }
}

my $bigaln = Bio::SimpleAlign->new;
while( my ($id,$seq) = each %matrix ) {
    $bigaln->add_seq(Bio::LocatableSeq->new(-id  => $id,
					    -seq => $seq));
}

my $out = Bio::AlignIO->new(-format => $oformat,
			    -file   => ">$outfile");
$out->write_aln($bigaln);


