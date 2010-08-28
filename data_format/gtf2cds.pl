#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;
my $prefix = '';
GetOptions(
	   'p|prefix:s' => \$prefix);
my $db = Bio::DB::Fasta->new(shift @ARGV);

if( $prefix && $prefix !~ /\.$/ ) {
    $prefix .= ".";
}

my %genes;
while(<>) {
    my @line = split(/\t/,$_);
    next unless $line[2] eq 'CDS';
    my %dat;
    for my $field ( qw(gene_id transcript_id exon_number) ) {
	if( $line[-1] =~ /$field\s+(\S+)/ ) {
	    $dat{$field} = $1;
	}
    }
    my $group = $dat{'transcript_id'};
    $group =~ s/[\";]//g;
    my $strand = 1;
    if($line[6] eq '-' ) {
#	warn("flip flopping @line\n");
	($line[3],$line[4]) = ($line[4],$line[3]);
	$strand = -1;
    }
    push @{$genes{$group}}, [$line[3],$line[4],$strand,$line[0],\%dat ];
}
my $outi = Bio::SeqIO->new(-format => 'fasta',
			   -file => '>'.$prefix.'introns.fa');
my $oute = Bio::SeqIO->new(-format => 'fasta',
			   -file => '>'.$prefix.'cds.fa');

for my $gene ( sort keys %genes ) {
  my ($cds,$lastexon);
  my $intron = 0;
  for my $exon ( map { $_->[0] }
		 sort { $a->[1] <=> $b->[1] }
		 map { [$_, $_->[2]*$_->[0]] }
		 @{$genes{$gene}} ) {
    $cds .= $db->seq($exon->[3], $exon->[0] => $exon->[1]);

    # this is the intron extraction part
    if ( $lastexon ) {
      my ($s,$e);
      # check what strand the exon is on
      if ( $exon->[2] < 0 ) {	# reverse strand exon 
	($s,$e) = ($lastexon->[1]-1,$exon->[0]+1);
      } else {			# fwd strand exon
	($s,$e) = ($lastexon->[1]+1, $exon->[0]-1);
      }
      # write out the exon with the header
      # >GENE.iINTRONNUM gene_id=GENE"
      $outi->write_seq(Bio::PrimarySeq->new
		       (-id => sprintf("%s.i%d",$gene,$intron++),
			-description => 
			sprintf("gene_id=%s",$exon->[4]->{'gene_id'}),
			-seq => $db->seq($exon->[3],$s => $e)));
    }
    $lastexon = $exon;
  }
  $oute->write_seq(Bio::PrimarySeq->new
		   (-id => $gene,
		    -seq => $cds));
}
