#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

=head1 USAGE

 four_gamete_test -i snpfile

Display this help
 four_gamete_test --help 

=head1 DESCRIPTION

 -l/--locus [fixed or gene] will specify whether to use fixed-size bins 
       or gene bins


=cut

use List::Util qw(shuffle sum);
use Getopt::Long;

my $snpfile = '/projects/coccidioides/variation/Cocci_integrated_SNPreport_6Aug09.txt';

my $locus = 'fixed';
my $Window_size = 1_000;
my $Group = 'CP'; # C.posa
my $Minimum_distance = 500; # so this is _bin_ distance so mult by
                            # Window_size to get bp distance
my $Num_loci = 500; # number of loci to consider
my $Min_locus_sites = 3;
GetOptions(
	   'l|locus:s'       => \$locus,
	   'i|input:s'       => \$snpfile,
	   'w|window:i'      => \$Window_size,
	   'm|min:i'         => \$Minimum_distance,
	   'n|num|numloci:i' => \$Num_loci,
	   'g|group:s'       => \$Group,
	   'h|help'    => sub { exec('perldoc',$0);
				exit(0);
			    },
	   );

open( my $fh => $snpfile) || die $!;
my $i = 0;
# make a hash out of the header
my %header = map { $_ => $i++ } split(/\s+/,<$fh>);

my %names = ( 'CP' => [grep { /^CP_/ } keys %header],
	      'CI' => [grep { /^CI_/ } keys %header]);
$names{ALL} = [ map { @{$names{$_}} } keys %names];
# keys will be first chromosomes and then we'll use an array to represent each
# this should all fit in memory, it is < 1M SNPs
my %chroms;

while(<$fh>) {
    chomp; # remove trailing '\n'
    my @row = split(/\t/,$_);
    
    my $chrom = $row[ $header{'supercontig'} ] ;
    next unless $row[ $header{'class'}] eq 'Intergenic';
    push @{$chroms{ $chrom }}, [@row];   
}

&calc_IA(\%chroms);

sub calc_IA {
    my $chroms = shift;
    my %loci;
    my %K;
    my $no_individuals = scalar @{$names{$Group}};
    for my $chrom ( sort { $a <=> $b } keys %$chroms ) {
	my $bin;
	my $dat = $chroms->{$chrom};
	for my $row ( @{$dat} ) {	
	    if( $locus eq 'gene') { 	
                # bin-by-gene
		next unless( defined ( $bin = $row->[$header{'gene'}] ) ); 
	    } elsif( $locus eq 'fixed' ) {
		# grab a sub-partition
                # bin-by-gene
		next unless( defined ( $bin = $row->[$header{'pos(bp)'}] ) ); 
		$bin /= $Window_size;
		$bin = int($bin);
	    } else {
		die "unknown locus type\n";
	    }

	    my %alleles;
	    my $total = 0;
	    for my $al ( grep { defined $_ && length($_) } 
			 map { $row->[ $header{$_} ] } 
			 @{$names{$Group}} ) {
		$alleles{$al}++;
		$total++;
	    }
	    next if( (scalar keys %alleles) <= 1 ); # skip monoallelic

	    for( my $i =0; $i < $no_individuals; $i++ ) {
		my $ind_i = $names{$Group}->[$i];
		for( my $j = $i+1; $j < $no_individuals; $j++ ) {
		    my $ind_j = $names{$Group}->[$j];
		    $K{"$ind_i,$ind_j"} ||= 0;
		    if( defined $row->[ $header{$ind_i} ] &&
			defined $row->[ $header{$ind_j} ] &&
			( $row->[ $header{$ind_j} ] ne 
			  $row->[ $header{$ind_i} ]) ){
			push @{$K{"$chrom.$bin"}}, "$ind_i,$ind_j";
		    }
		}
	    }
	    my $sum;
	    for my $r ( keys %alleles ) {
		my $pij = ($alleles{$r} /= $total);		
		$sum += ($pij ** 2);
		#warn("freq is $alleles{$r}\n");
	    }
	    push @{$loci{$chrom.".".$bin}}, $sum; 
	}
    }
    my @locus_names = keys %loci;
    warn(sprintf("there are %d loci\n",scalar @locus_names));
    my $K_mean;			# mean difference between individuals
    my $Var_E;			# expected variance

# this is the part where we need to subsample

    
    my $count;
    my %used;
    my %final_K;
    for my $bin ( shuffle @locus_names ) {
	if( &too_close($bin, \%used)) {
	    # warn("$bin is too close to previously seen bins\n");
	    next;
	}
	$used{$bin}++;
	next unless( scalar @{$loci{$bin}} >= $Min_locus_sites);
	for my $pair ( @{$K{$bin}} ) {
	    $final_K{$pair}++;
	}
	my $hj = 1 - sum(@{$loci{$bin}});
	$K_mean += $hj;
	$Var_E  += ($hj * ( 1 - $hj));
	last if $count++ > $Num_loci;
    }
    my $stats = Statistics::Descriptive::Full->new;
    $stats->add_data(values %final_K);
    
    warn("count is $count loci used\n");
    printf "V_O = %.2f\n", $stats->variance();
    printf "V_E = %.2f\n", $Var_E;
    printf "I_A = %.2f\n", $stats->variance()/$Var_E -1;

}
sub too_close {
    my ($locus, $used_loci) = @_;
    my ($cur_chr,$cur_pos) = split(/\./,$locus);
    for my $l ( keys %$used_loci ) {
	my ($chr,$pos) = split(/\./,$l);
	if( $cur_chr == $chr &&
	    (abs($cur_pos - $pos) < $Minimum_distance)) {
	    #warn("  $locus is too close to $l\n");
	    return 1;
	}
    }
    return 0;
}
	       

sub bootstrap {
    # shuffle the allele assignments
}
