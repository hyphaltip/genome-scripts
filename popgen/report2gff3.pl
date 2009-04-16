#!/usr/bin/perl -w
use List::Util qw(sum);
use strict;
use Getopt::Long;
my $debug = 0;
my $src = 'BroadSNP_13';
my $Datafile = '/projects/coccidiodies/variation/SnpReport_Cspp-3.txt';

GetOptions('v|verbose|debug!' => \$debug,
	   's|src:s'          => \$src,
	   'd|data|i|in:s'    => \$Datafile,
	   );

open( my $fh => $Datafile ) || die "cannot open $Datafile: $!";
my $i = 0;
chomp(my $hdr = <$fh>);
$hdr =~ s/^\#//;
my %header = map { $_ => $i++ } split(/\t/,$hdr);

# A speedup to pre-index the by-species allele columns
my (@allele_cols,@sp_names, %sp_count, @HdrCols);

print "##gff-version 3\n";

for my $col ( grep { /^a_(\S+)/ } keys %header ) {
    unless ( $col =~ /^a_(C[IP])/ ) {
	warn("unknown allele column '$col'\n");
    }
    $sp_count{$1}++;
    $sp_count{'ALL'}++; # for total number of sample count, when data is not partitioned by species
    push @allele_cols, [ $1, $header{$col} ];
}

@sp_names = sort grep { !/^ALL$/ } keys %sp_count;
while(<$fh>) {
    chomp($_); # remove the trailing '\n'
    my @row = split(/\t/,$_); # split into an array
    my ($snpid,$lg,$start, $refallele,
	$alleles,$genes) = 
	    map { $row[ $header{$_} ] } qw(snpid linkageGroup pos 
					   refallele alleles genes);
	    
    $lg = 'cimm_RS_2.'.$lg;

    # if we did re-sampling/removal of strains, we'd do it here
    my $allele_freq = {};     
    for my $col ( @allele_cols ) {
	my $allele = $row[ $col->[1] ];
	next unless length($allele);
	$allele_freq->{$col->[0]}->{$allele}++;	
	$allele_freq->{'ALL'}->{$allele}++;  # for the all (Group_by_species = false) combined data		
    }    
    my @acount;
    for my $sp ( @sp_names ) {
	my @allele_info;
	# total count of alleles in this sp
	my $total = sum ( values %{$allele_freq->{$sp}} );
	for my $allele ( split('/',$alleles) ) {
	    my $n = exists $allele_freq->{$sp}->{$allele} ? $allele_freq->{$sp}->{$allele} : 0;
	    my $num; # format so it is 0,1, or a 0.50 ... etc
	    if( $n == 0 ) {
		$num = 0;
	    } elsif( $n == $total ) {
		$num = 1;
	    } else {
		$num =sprintf("%.3f", $n / $total);
	    }
	    push @allele_info, sprintf("%s+%s+%d",$allele,$num,$n);
	}
	push @acount, sprintf("%s:%s",$sp,join(" ",@allele_info));
    }    
    print join("\t",
	       $lg, $src, 'snp', $start, $start, '.', '.','.',
	       sprintf("ID=SNP:%s;Name=%s;Alleles=%s;Refallele=%s;Acounts=%s;Gene=%s",
		       $snpid,$snpid,$alleles,$refallele,
		       join(",",@acount),$genes)
	       ),"\n";
    last if $debug;
}
