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

GetOptions('d|dir:s'   => \$dir,
	   'if|informat:s' => \$iformat,
	   'of|outformat:s' => \$oformat,
	   'ext:s'     => \$ext,
	   'o|out:s'   => \$outfile,
	   );

die("need a dir") unless $dir && -d $dir;

opendir(DIR, $dir) || die"$dir: $!";

my (%matrix);

for my $file (sort readdir(DIR) ) {
    next if $file eq $outfile;
    next unless ($file =~ /(\S+)\.\Q$ext\E$/);
    my $in = Bio::AlignIO->new(-format => $iformat,
			       -displayname_flat => 1,
			       -file   => "$dir/$file"); 
     my %seen;
    if( my $aln = $in->next_aln ) {
	for my $seq ( $aln->each_seq ) {
	    if (length($seq->seq) != $aln->length ) {
		warn($seq->id, " problem with length in ", $file,"\n");
	    }
	    $matrix{$seq->id} .= $seq->seq;
	}
    }
}
my $bigaln = Bio::SimpleAlign->new;
while( my ($id,$seq) = each %matrix ) {
    warn("$id ", length($seq),"\n");
    $bigaln->add_seq(Bio::LocatableSeq->new(-id  => $id,
					    -seq => $seq));
}

my $out = Bio::AlignIO->new(-format => $oformat, -show_symbols => 0, -show_endblock=> 0,
			    -file   => ">$outfile");
$out->write_aln($bigaln);

