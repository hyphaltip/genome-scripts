#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Env qw(USER HOME);
use Getopt::Long;
use File::Spec;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use Data::Dumper;

my $debug;
my $dir = 'promoter_alignments';
my @sp  = qw(Ncra Ntet Ndis);
my $genome = 'Smac';
GetOptions(
	   'i|in|dir:s' => \$dir,
	   'g|genome:s' => \$genome,
	   'v|verbose!'=> \$debug,
	   );

my $stats = Bio::Align::DNAStatistics->new;

opendir(DIR, $dir) || die $!;

for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.fasta$/;
    my $stem = $1;
    my $alndir = File::Spec->catfile($dir,"$stem.d");
    opendir(ALNDIR, $alndir) || die $!;
    open(my $fh => ">$dir/$stem.stats.dat") || die $!;
    warn("trying out $stem.stats.dat\n");
    print $fh join("\t", qw(FEATURE ALN_LENGTH SEQ_LEN ALL_ID NC_JCUC PW_ID), 
		   map { uc($_ ."_D") } @sp),"\n";
    for my $fd ( readdir(ALNDIR) ) {
	next unless ( $fd =~ /(\S+)\.fas$/);
	my $fstem = $1;
	my $in = File::Spec->catfile($alndir,$fd);
	my $input = Bio::AlignIO->new(-alphabet => 'dna',
				      -format => 'fasta',
				      -file   => $in);
	my $aln = $input->next_aln;	
	next unless defined $aln;
	my $skip;
	for my $s ( $aln->each_seq ) {
	    $skip = 1 if ( $s->end == 0 );
	}
	next if $skip;
	my $length = $aln->length;
	my $dist;
	eval {
	    $dist = $stats->distance(-align => $aln,
				     -method => 'k80');	
	};
	if( $@ ) {
	    warn("$@ for $in\n");
	}
	next unless( $dist && defined $dist);       
	
	my $new_pw_aln = Bio::SimpleAlign->new;
	my $ncra = $aln->get_seq_by_id('Ncra');
	my $smac = $aln->get_seq_by_id('Smac');
	my $smac_seq = $smac->seq;
	$smac_seq =~ s/\-//g;
	my $seq_len = length($smac_seq);
	my ($pw_id,$pw_JC) = (0,0);
	if( $ncra->end > 0 ) {
	    $new_pw_aln->add_seq($ncra);
	    $new_pw_aln->add_seq($smac);
	    $new_pw_aln->remove_gaps(1);

	    $pw_id = $new_pw_aln->overall_percentage_identity;
	    if( $new_pw_aln->length > 0 && $pw_id > 0 ) {
		my $uncorr_dist = $stats->distance(-align => $aln,
						   -method => 'uncorrected');
		$pw_JC = 100 * (1 - $uncorr_dist->get_entry($genome,'Ncra'));
	    }
	}
	print $fh join("\t", $fstem, $length,$seq_len,
		       sprintf("%.1f",$aln->average_percentage_identity),
		       sprintf("%.2f",$pw_JC),sprintf("%.1f",$pw_id),
		       map { $dist->get_entry($genome,$_) } @sp),"\n";
		   
	last if $debug;
    }
    last if $debug;
}
