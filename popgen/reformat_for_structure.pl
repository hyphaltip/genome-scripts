#!/usr/bin/perl -w
my $k = 0;
my %base2num = map { $_ => $k++ } qw(N A C G T);
use strict;
use List::Util qw(sum);
use Getopt::Long;
my $debug = 0;
my $Datafile = '/projects/coccidiodies/variation/SnpReport_Cspp-3.txt.bz2';
my $all_individuals = 0;
my $population = 'CP';
my $outdir  = 'structure';
my $combine = 0;
GetOptions('v|verbose|debug!' => \$debug,
	   'd|data|i|in:s'    => \$Datafile,
	   'p|population:s'   => \$population,
	   'a|all|allinds!'   => \$all_individuals,
	   'o|out|output:s'   => \$outdir,
	   'c|combine!'       => \$combine,
	   );
mkdir($outdir) unless -d $outdir;
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
	    push @alleles, [$col, $base2num{$allele}];
	    $al{$base2num{$allele}}++;
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
    #$outgroup{$lg}->{$snpid} = $base2num{$refallele};
    last if $debug && $i++ > 1000;
}

if( $combine ) {
    open(my $rptfh => ">$outdir/SnpReport_Cspp-3.structure.tab") || die $!;
    my @chroms = sort keys %data;
    my (@distances,@all_snps,@all_rows);
    for my $chrom ( @chroms ) {	
	my @snps = sort { $a->[1] <=> $b->[1] } map { [$_,(split('_',$_))[-1]] } keys %{$data{$chrom}};
	my $ct = scalar @snps;
	
	# print inter-marker distance
	my $lastsnp;
	for my $snp ( map { $_->[1] } @snps ) {
	    unless( $lastsnp ) {
		push @distances, -1;
	    } else {
		push @distances, $snp - $lastsnp;
	    }
	    $lastsnp = $snp;
	}
	push @all_snps, map { $_->[0] } @snps;
	my $i = 0;
	for my $ind ( @allele_cols ) {
	    push @{$all_rows[$i++]}, 
	    map {  $data{$chrom}->{$_->[0]}->{$ind} } @snps;
	}
    }
    print $rptfh join(" ", @all_snps),"\n";
    print $rptfh join(" ", @distances), "\n";
    my $i = 0;
    for my $ind ( @allele_cols ) {
	my $ind_print = $ind;
	$ind_print =~ s/^a_//;	
	print $rptfh join(" ",$ind_print, @{$all_rows[$i++]}),"\n";       
    }
    warn( (scalar @all_snps)," snps\n");
} else {

    for my $chrom ( sort keys %data ) {
	open(my $rptfh => ">$outdir/SnpReport_Cspp-3.$chrom\_structure.tab") || die $!;

	my @snps = sort { $a->[1] <=> $b->[1] } map { [$_,(split('_',$_))[-1]] } keys %{$data{$chrom}};
	@snps = splice(@snps,0,100) if $debug;
	my $ct = scalar @snps;    
	print $rptfh join(" ", map { $_->[0] } @snps),"\n";

	# print inter-marker distance
	my $lastsnp;
	my @distances;
	for my $snp ( map { $_->[1] } @snps ) {
	    unless( $lastsnp ) {
		push @distances, -1;
	    } else {
		push @distances, $snp - $lastsnp;
	    }
	    $lastsnp = $snp;
	}
	print $rptfh join(" ", @distances), "\n";

	for my $ind ( @allele_cols ) {
	    my $ind_print = $ind;
	    $ind_print =~ s/^a_//;
	    print $rptfh join(" ",$ind_print, map {  $data{$chrom}->{$_->[0]}->{$ind} } @snps), "\n";
	}
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
