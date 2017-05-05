#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
my $noref = 0;
my $control_strain; # expect this to be same as ref genome
my $ref_strain = 'REF';
my $invariants = 1; # include sites which are fixed in all strains?
my $ofilename;
GetOptions('noref!' => \$noref,
           'refstrain:s' => \$ref_strain,
           'I|invariant!' => \$invariants,
	   'o|ofile:s'   => \$ofilename,
	   'c|control:s' => \$control_strain);

if( $control_strain && $ref_strain eq 'REF' ) {
 $ref_strain = $control_strain;
}
if (! $ofilename && @ARGV) {
 $ofilename = $ARGV[0];
 $ofilename =~ s/\.tab//;
} elsif ( ! $ofilename ) {
 $ofilename = 'vcf2fastaseq';
}
my $head = <>;
my ($chrom,$pos,$ref,@strains) = split(/\s+/,$head);

my $i = 0;
my %strain2col = map { $_ => $i++ } @strains;
my $total = scalar @strains;
my %seqs;
my $ctl_col = $control_strain && exists $strain2col{$control_strain} ? $strain2col{$control_strain} : -1;

while(<>) {
    my ($chr,$p,$ref_allele, @alleles) = split;
    my %count = ( $ref_allele => 1 );
    my @sanitized;
    for my $allele ( @alleles ) {
	my ($al1,$al2) = split(/\//,$allele);
        
	if( $al2 && $al2 ne $al1 ) {
	    # this is diploid
	    # do something different?
	    warn " cannot process diploid tab files for now ($allele)\n";
	} 
	$al1 =~ s/\./-/;
	$count{$al1}++;
	push @sanitized, $al1;
    }

    if( exists $count{'-'} && 
	(my $per = $count{'-'} / $total) > 0.9 ) {
	warn(sprintf("skipping $chr,$p as %.2f%% are uncalled\n", 
		     $per*100));
	next;
    }
    if( $ctl_col >= 0 && $ref_allele ne $sanitized[$ctl_col] ) {
	warn("Ref @ $chr:$p = $ref_allele but reseq strain allele is ",$sanitized[$ctl_col]," replacing\n");
        $seqs{$ref_strain} .= $sanitized[$ctl_col] unless $noref;
    } else {
     $seqs{$ref_strain} .= $ref_allele unless $noref;
    }
    my $n = 0;
    my %af;
    for my $al ( @sanitized ) {
	$af{$al}++ if ( $n++ != $ctl_col );
    }

    next if ( ! $invariants && (scalar keys %af) == 1 ); # skip when site is invariant unless you ask for invariants (default)
    $n = 0;
    for my $al1 ( @sanitized ) {
	if( $n != $ctl_col ) {
	 $seqs{$strains[$n]} .= $al1;	
        }
        $n++;
    }
}
if( $ofilename !~ /\.(fas|fasta|seq|nt|dna)/i ) {
 $ofilename .= ".fas";
}
my $out = Bio::SeqIO->new(-format => 'fasta',-file => ">$ofilename");
for my $strain ( keys %seqs ) {
    $out->write_seq(Bio::Seq->new(-id => $strain,
				  -seq => $seqs{$strain}));
}
