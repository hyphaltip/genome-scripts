#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;
my $db = Bio::DB::Fasta->new(shift @ARGV);
my $prefix = '';
GetOptions(
	   'p|prefix:s' => \$prefix);

if( $prefix && $prefix !~ /\.$/ ) {
    $prefix .= ".";
}

my (%genes, @order);
while(<>) {
    my @line = split(/\t/,$_);
    next if $line[2] eq 'mRNA';
    my (undef,$group) = split(/\s+/,$line[-1]);
    $group =~ s/\"//g;
    my $strand = 1;
    if($line[6] eq '-' ) {
#	warn("flip flopping @line\n");
	($line[3],$line[4]) = ($line[4],$line[3]);
	$strand = -1;
    }
    push @order, $group unless exists $genes{$group};
    push @{$genes{$group}}, [$line[3],$line[4],$strand,$line[0] ];
}
my $outi = Bio::SeqIO->new(-format => 'fasta',
			   -file => '>'.$prefix.'introns.fa');
my $oute = Bio::SeqIO->new(-format => 'fasta',
			   -file => '>'.$prefix.'cds.fa');
my $outp = Bio::SeqIO->new(-format => 'fasta',
		           -file => '>'.$prefix.'pep.fa');
for my $gene ( @order ) {
    my ($cds,$lastexon);
    my $intron = 0;
    for my $exon ( map { $_->[0] }
		    sort { $a->[1] <=> $b->[1] }
		    map { [$_, $_->[2]*$_->[0]] }
		    @{$genes{$gene}} ) {
	$cds .= $db->seq($exon->[3], $exon->[0] => $exon->[1]);
	if( $lastexon ) {
	    my ($s,$e);
	    if( $exon->[2] < 0 ) {
		($s,$e) = ($lastexon->[1]-1,$exon->[0]+1);
	    } else {
		($s,$e) = ($lastexon->[1]+1, $exon->[0]-1);
		
	    }
	    $outi->write_seq(Bio::PrimarySeq->new
			     (-id => sprintf("%s.i%d",$gene,$intron++),
			      -seq => $db->seq($exon->[3],
					       $s => $e)));
	}
	$lastexon = $exon;
    } 
    my $s = Bio::PrimarySeq->new(-id => $gene, -seq => $cds);
    $oute->write_seq($s);
    $outp->write_seq($s->translate);
}
