#!/usr/bin/perl -w
use strict;

my $Datafile = shift || 'SnpReport_Cspp-3.txt';
open(my $ci_snp  => ">plot/snp_ci_specific.tab") || die $!;
open(my $cp_snp  => ">plot/snp_cp_specific.tab") || die $!;
open(my $both_snp  => ">plot/snp_both_poly.tab") || die $!;
open(my $all_snp => ">plot/snp_all.tab") || die $!;
open(my $div_snp => ">plot/snp_divergences.tab") || die $!;

for my $fh ( $ci_snp, $cp_snp, $all_snp, $div_snp,$both_snp ) {
    print $fh '#', join("\t",  qw(src src_start ID)),"\n";
}


open(my $in => $Datafile ) || die $!;
my $i = 0;
chomp(my $hdr = <$in>);
$hdr =~ s/^\#//;
my @header = split(/\t/,$hdr);
my %header = map { $_ => $i++ } @header; 

# A speedup to pre-index the by-species allele columns
my (@allele_cols,@sp_names, %sp_count);

for my $col ( grep { /^a_(\S+)/ } keys %header ) {
    unless ( $col =~ /^a_(C[IP])/ ) {
	warn("unknown allele column '$col'\n");
    }
    $sp_count{$1}++;
    $sp_count{'ALL'}++; # for total number of sample count, when data is not partitioned by species
    push @allele_cols, [ $1, $header{$col} ];
}
@sp_names = sort grep { !/^ALL$/ } keys %sp_count;

while(<$in>) {
    chomp($_); # remove the trailing '\n'
    my @row = split(/\t/,$_); # split into an array
    my @data = map { $row[ $header{$_} ] } qw(linkageGroup pos snpid);
    $data[0] = 'cimm_2.'.$data[0];
    print $all_snp join("\t", @data),"\n";
    # if we did re-sampling/removal of strains, we'd do it here
    my %allele_freq;
    my %af;
    for my $col ( @allele_cols ) {
	my $allele = $row[ $col->[1] ];
	next unless length($allele);
	$allele_freq{$col->[0]}->{$allele}++;	
	$af{$col->[0]}->{$header[$col->[1]]} = $allele;
    }

    if( keys %allele_freq == 1  ) {  # only 1 species has a SNP
	if( exists $allele_freq{'CP'} ) { # this is a divergence
	    if( keys %{$allele_freq{'CP'}} > 1 ) {
		# this is a CP polymorphism
		print $cp_snp join("\t", @data), "\n";
	    } else {
		# this will never happen
		print $div_snp join("\t", @data), "\n";
	    }
	} else { # this is a Ci polymorphism
	    print $ci_snp join("\t", @data), "\n";
	}
    } else {
	# if there are polymorphisms in both species this 
	# is a CI and CP SNP
	if( keys %{$allele_freq{'CP'}} > 1 && 
	    keys %{$allele_freq{'CI'}} > 1 ) {
	    print $both_snp join("\t", @data),"\n";
#				 join(",", 
#				map { sprintf("%s=%s",$_,$af{'CI'}->{$_})}
#				      keys %{$af{'CI'}}),
#				map { sprintf("%s=%s",$_,$af{'CP'}->{$_})}
#				 keys %{$af{'CP'}}),"\n";
	} elsif( keys %{$allele_freq{'CI'}} > 1 )  {
	   print $ci_snp join("\t", @data), "\n"; 
	} elsif( keys %{$allele_freq{'CP'}} > 1 ) {
	   print $cp_snp join("\t", @data), "\n"; 
       } else {
	   print $div_snp join("\t", @data), "\n";
       }
	
    }
}
