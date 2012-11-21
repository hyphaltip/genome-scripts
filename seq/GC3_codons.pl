#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(sum);
use Bio::Tools::CodonTable;

my $pergene;

GetOptions('g|pergene!' => \$pergene);

if( $pergene ) {
 print join("\t",qw(SEQ GC G C)), "\n";
}
my $infile = shift;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => $infile);

my $codontable = Bio::Tools::CodonTable->new();
my %aa;

my @bases = qw(A C G T);
for my $firstbase ( @bases ) {
  for my $secondbase ( @bases ) {
    for my $thirdbase ( @bases ) {
      my $codon = $firstbase . $secondbase . $thirdbase;
      push @{$aa{$codontable->translate($codon)}}, $codon;
    }
  }
}
my %degenerate;
while ( my ($aa,$codons) = each %aa ) {
  if ( @$codons >= 4 ) { # 4fold or 6fold degenerate
    for my $codon ( @$codons)  {
      $degenerate{$codon}++;
    }
  }
}
my (%gc,%nongc);
while (my $seq = $in->next_seq ) {
  my $str = $seq->seq;
  if ( ($seq->length % 3) != 0 ) {
    warn($seq->id, " coding sequence is not multiple of 3 ", $seq->length, "\n");
    next;
  } else {
#    warn($seq->id, " ok\n");
  }
  my (%gene_count,%nongene_count);
  while ( $str )  {
   my $codon = substr($str,0,3,'');
   my $thirdbase = uc substr($codon,2,1);
   if ( $degenerate{$codon} ) {
     $gc{$thirdbase}++;
     $gene_count{$thirdbase}++;
   } else {
     $nongc{$thirdbase}++;
     $nongene_count{$thirdbase}++;
   }
 }
  next if ! keys %gene_count;
  if ( $pergene ) {
	my $totalbp = sum(map { $gene_count{$_} || 0} @bases);	
	my $g = $gene_count{'G'} || 0;
	my $c = $gene_count{'C'} || 0;
	my $non_totalbp = sum(map { $nongene_count{$_} || 0} @bases);	
	my $non_g = $nongene_count{'G'} || 0;
	my $non_c = $nongene_count{'C'} || 0;
	
    printf "%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",$seq->id, ($g+$c) / $totalbp, $g/$totalbp, $c/$totalbp, 
	($non_g + $non_c) / $non_totalbp,
	$non_g / $non_totalbp, $non_c/$non_totalbp;
  }
}

my $gc = ($gc{'G'} || 0)  + ($gc{'C'} || 0);
my $g = $gc{'G'};
my $c = $gc{'C'};
my $total = sum ( map { $gc{$_} } @bases) ;
print join("\t", qw(CODON_TYPE GC G C)), "\n";
printf "%s\t%.3f\t%.3f\t%.3f\n",'4+6 fold DG',($g+$c) / $total, $g/$total, $c/$total;

$gc = ($nongc{'G'} || 0)  + ($nongc{'C'} || 0);
$g = $nongc{'G'};
$c = $nongc{'C'};
$total = sum ( map { $nongc{$_} } @bases);

printf "%s\t%.3f\t%.3f\t%.3f\n",'non 4+6 fold DG', ($g+$c) / $total, $g/$total, $c/$total;
