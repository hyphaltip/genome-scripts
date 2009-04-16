#!/usr/bin/perl -w

use strict;
use List::Util qw(sum);
use Getopt::Long;
my $debug = 0;
my $Datafile = '/projects/coccidiodies/variation/SnpReport_Cspp-3.txt.bz2';
my $all_individuals = 0;
my $population = 'CP';
my $rsq = 'rsq';
GetOptions('v|verbose|debug!' => \$debug,
	   'd|data|i|in:s'    => \$Datafile,
	   'p|population:s'   => \$population,
	   'a|all|allinds!'   => \$all_individuals,
	   'o|output:s'       => \$rsq,
	   );
mkdir($rsq) unless -d $rsq;
open( my $fh => "bzcat $Datafile |" ) || die "cannot open $Datafile: $!";
my $i = 0;
chomp(my $hdr = <$fh>);
$hdr =~ s/^\#//;
my %header = map { $_ => $i++ } split(/\t/,$hdr);

# A speedup to pre-index the by-species allele columns
my (@allele_cols,@sp_names, %sp_count, @HdrCols);

for my $col ( grep { /^a_(\S+)/ } keys %header ) {
    unless ( $col =~ /^a_(C[IP])/ ) {
	debug("unknown allele column '$col'\n");
    }
    my $sp = $1;
    next unless $sp eq $population;
    push @allele_cols, $col;
}
$i =0;
my (%outgroup,%data);
while(<$fh>) {
    chomp($_); # remove the trailing '\n'
    my @row = split(/\t/,$_); # split into an array
    
    my ($snpid,$lg,$refallele) = map { $row[ $header{$_} ] } qw(snpid linkageGroup refallele);
    $lg = 'cimm_RS_2.'.$lg;
    
    my (@alleles,%al,%inds);
    for my $col ( @allele_cols ) {
	my $allele = $row[ $header{$col} ];
	if( $allele && length($allele) ) {
	    push @alleles, [$col,$allele];
	    $al{$allele}++;
	    $inds{$col}++;
	}
    } 
    # insure this is bi-allelic
    next unless( (keys %al) == 2 );
    next if( $all_individuals && (keys %inds) != @allele_cols );
    for my $al ( @alleles ) {
	my ($ind, $allele) = @$al;
	$data{$lg}->{$snpid}->{$ind} = $allele;
    }
    $refallele = 'N' if $refallele eq '-' || ! $refallele;
    $outgroup{$lg}->{$snpid} = $refallele;
    last if $debug && $i++ > 1000;
}

for my $chrom ( sort keys %data ) {
    open(my $rptfh => ">$rsq/SnpReport_Cspp-3.$chrom\_summary.txt") || die $!;
    
    my @snps = keys %{$data{$chrom}};
    my $ct = scalar @snps;
    printf $rptfh "%d %d\n",scalar @allele_cols, $ct;
    print $rptfh join(" ", map { int((split('_',$_))[-1]) } @snps),"\n";
    print $rptfh join(" ", map { $outgroup{$chrom}->{$_} } @snps),"\n";
    for my $ind ( @allele_cols ) {
	print $rptfh join(" ",$ind, map { $data{$chrom}->{$_}->{$ind} || 'N' } @snps),
	"\n";
    }
}

sub intersection {
    my @lsts = @_;
    my $expected = scalar @lsts;
    my %lst;
    $lst{$_}++ for ( map { @$_ } @lsts );
    return sort grep { $lst{$_} == $expected } keys %lst;
}


sub debug {
    my $line = shift;
    warn($line) if $debug;
}
