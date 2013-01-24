#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use List::Util qw(sum);
use Bio::Tools::CodonTable;

my $pergene;
my $min_length = 3;
GetOptions('g|pergene!' => \$pergene,
	   'l|length:i' => \$min_length);

# print out the header about the sequence info
if( $pergene ) {
 print join("\t",qw(SEQ GC3 G3 C3 GC_ALL G_ALL C_ALL)), "\n";
}
my $infile = shift;

# sequence parser
my $in = Bio::SeqIO->new(-format => 'fasta', -file => $infile);

# us this to translate a codon into its amino acid
# with this we can map multiple codons to the same amino acid, and thus determine which ones 
# have synonymous codons
my $codontable = Bio::Tools::CodonTable->new();
my %aa;

# this just enumerates all the codons (all combinations of 3 letters of DNA
my @bases = qw(A C G T);
for my $firstbase ( @bases ) {
  for my $secondbase ( @bases ) {
    for my $thirdbase ( @bases ) {
      my $codon = $firstbase . $secondbase . $thirdbase;
      push @{$aa{$codontable->translate($codon)}}, $codon;
    }
  }
}

# this will make the hash %degerate which will have a 1 if the codon is a 4 or 6 fold degenerate codon
# and will have no value if it is not, so we can test last, given a codon, is it 4/6-fold degenerate
my %degenerate;
while ( my ($aa,$codons) = each %aa ) {
  if ( @$codons >= 4 ) { # 4fold or 6fold degenerate
    for my $codon ( @$codons)  {
      $degenerate{$codon}++;
    }
  }
}
my %gc = map { $_ => 0 } @bases;
my %nongc = map { $_ => 0} @bases;

while (my $seq = $in->next_seq ) {
    next if $seq->length < $min_length;
    my $str = $seq->seq;
    # warn if this sequence is not divisible perfectly by 3, so it is a partial gene prediction
    if ( ($seq->length % 3) != 0 ) {
	warn($seq->id, " coding sequence is not multiple of 3 ", $seq->length, "\n");
	next;
    } else {
#    warn($seq->id, " ok\n");
    }
    # initialize variables for the nucleotide counts
    
    my %deg_count = map { $_ => 0 } @bases;
    my %nondeg_count = map { $_ => 0 } @bases; # degenerate codons count, and non-degenerate codons count
    
    while ( $str )  {
	my $codon = uc substr($str,0,3,''); # pull first 3 letters off the string (and remove them from the string
	my $thirdbase = substr($codon,2,1); # 3rd letter is position 2 in the string
	next if $thirdbase eq 'N';
	if ( $degenerate{$codon} ) {
	    $gc{$thirdbase}++;
	    $deg_count{$thirdbase}++;
	} else {
	    $nongc{$thirdbase}++;
	    $nondeg_count{$thirdbase}++;
	}
    }
    next if ! keys %deg_count;
    if ( $pergene ) {
	my $totalbp = sum(values %deg_count );
	my $g = $deg_count{'G'};
	my $c = $deg_count{'C'};
	my $non_totalbp = sum(values %nondeg_count);
	my $non_g = $nondeg_count{'G'};
	my $non_c = $nondeg_count{'C'};

	print join ("\t", $seq->id, sprintf("%.3f",($g+$c) / $totalbp), 
		    sprintf("%.3f",$g/$totalbp), sprintf("%.3f",$c/$totalbp), 
		    sprintf("%.3f",($non_g + $non_c) / $non_totalbp),
		    sprintf("%.3f",$non_g / $non_totalbp), 
		    sprintf("%.3f",$non_c/$non_totalbp)), "\n";
    }
}

my $gc = $gc{'G'} + $gc{'C'};
my $g = $gc{'G'};
my $c = $gc{'C'};
my $total = sum ( values %gc);
warn join("\t", qw(CODON_TYPE GC G C)), "\n";
warn sprintf "%s\t%.3f\t%.3f\t%.3f\n",'4+6 fold DG',($g+$c) / $total, $g/$total, $c/$total;

my $nongc = $nongc{'G'} + $nongc{'C'};
my $nong = $nongc{'G'};
my $nonc = $nongc{'C'};
$total = sum ( values %nongc);

warn sprintf "%s\t%.3f\t%.3f\t%.3f\n",'non 4+6 fold DG', ($g+$c) / $total, $g/$total, $c/$total;
