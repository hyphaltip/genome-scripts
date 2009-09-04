#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(sum shuffle);

# see http://guanine.evolbio.mpg.de/lian/lianDoc/

my $snpfile = '/projects/coccidioides/variation/Cocci_integrated_SNPreport_6Aug09.txt';
my $Num_loci = 20;
my $Bin_size = 20_000;
my $Group = 'CP';
GetOptions(
	   'i|input:s'       => \$snpfile,
	   'b|bin:i'         => \$Bin_size,
	   'n|num|numloci:i' => \$Num_loci,
	   'g|Group:s'       => \$Group,
	   'h|help'    => sub { exec('perldoc',$0);
				exit(0);
			    },
	   );

open( my $fh => $snpfile) || die $!;
my $i = 0;
# make a hash out of the header
my %header = map { $_ => $i++ } split(/\s+/,<$fh>);
#warn("header has ", scalar keys %header, " cols\n");
my %names = ( 'CP' => [grep { ! /^CP_Silveira_ILLUMINA/ } grep { /^CP_/ } keys %header],
	      'CI' => [grep { /^CI_/ } keys %header]);
$names{ALL} = [ map { @{$names{$_}} } keys %names];
# keys will be first chromosomes and then we'll use an array to represent each
# this should all fit in memory, it is < 1M SNPs

open(my $out => ">$Group.dat") || die $!;
my @names_lst = @{$names{$Group}};
print $out join("\t", qw(SNP CHROM POS), @names_lst), "\n";
my %bins;
while(<$fh>) {
    my @row = split(/\t/,$_);
    chomp($row[-1]);
    
    my ($snpid,$chrom,$pos) = map { $row[ $header{$_} ] } qw(snpid supercontig pos(bp));
    my %alleles;
    my ($skip,$total) = (0,0);
    my @alleles_order;
    for my $nm ( @names_lst ) {
#grep { defined $_ && length($_) } 
	#map { $row[ $header{$_} ] } 
	my $al =$row[ $header{$nm} ];
	if( ! defined $al || ! length($al) ) {
	    $skip = 1;
	    #warn("missing an allele for $nm ($snpid)\n");
	    last;
	}
	$alleles{$al}++;
	push @alleles_order,$al;
	$total++;
    } 

    next if( $skip );
    next if (scalar keys %alleles <= 1 ); # skip monoallelic?
                                                    # SEE HERE 1
    if( ! defined $pos || ! defined $snpid ) {
	warn("line is empty $_");
    }
    push @{$bins{$chrom}->[int($pos / $Bin_size)]}, [$pos,$total,\@alleles_order,$snpid];
    print $out join("\t",$snpid, $chrom, $pos, @alleles_order), "\n";
}
close($out);
my %final = map { $_ => [] } @names_lst;


for my $chrom ( keys %bins ) {
    my $r = 0;
    for my $bin ( shuffle @{$bins{$chrom}} ) {
	last if @{$final{$names_lst[0]}} >= $Num_loci;
	next if $r++ % 2 == 1; # alternate
	my ($first,@loci) = shuffle( @$bin );
	# just take the 1st locus, it doesn't matter since it is shuffled
	# basically we are just flipping from vertical to horizonal
	# so we have to unwind
	for( my $j=0;$j < scalar @names_lst; $j++ ) {
	                         # GO TO SEE HERE 1 for the order
	    my $allele = $first->[2]->[$j];
	    next unless defined $allele;
	    push @{$final{$names_lst[$j]}}, $allele;
	} 
	
    }
}
warn("total loci is ", scalar @{$final{$names_lst[0]}}, "\n");
open($out => ">$Group.lian.in") || die $!;
for my $nm ( @names_lst ) {
    print $out join(" ", $nm, @{$final{$nm}}), "\n";
}
