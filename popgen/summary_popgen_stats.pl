#!/usr/bin/perl -w
# $Id: summary_popgen_stats.pl 250 2008-12-12 22:20:20Z stajich $
# 
use strict;
use Getopt::Long;
use Bio::PopGen::Statistics;
use List::Util qw(sum);

=head1 NAME
 summary_popgen_stats - script for calculating summary statistics for Coccidioides population data in Broad tab format 

=head1 USAGE

    [-m tajima_D] [-d DATAFILE] [-s l] [-w 1000] [--groupbyspecies (default) OR --nogroupbyspecies]

Options:
 -m or --method METHOD : tajima_D, fu_li_D, fu_li_F, pi, theta
 -d or --data DATAFILE : provide a different datafile than default
 -s or --strategy [w l]: by 'w'indow or 'l'ocus calculation, default is 'l'
 -w or --window SIZE   : for 'w' calculations, provide a different windowsize than default 
 -g or --groupbyspecies: boolean, process C.immitis and C.posas separately 
 --nog or 
     --nogroupbyspecies: prefix with 'no' to use all of them the 

=head1 AUTHORS

Jason Stajich, jason_stajich _AT_ berkeley.edu
 http://fungalgenomes.org/

=cut


# Hardcode for Tom S's file from the cocci project
my $genefile = 'plot/ci_gene_summary.tab';
my %genes;
{
    open(my $fh => $genefile) || die $!;
    my $header = <$fh>;
    $header =~ s/^\#//;
    my $i =0;
    my %header = map { $_ => $i++ } split(/\t/,$header);
    while(<$fh>) {
	my @row = split;
	$genes{$row[ $header{'gene'} ]} = [ $row[ $header{'chrom_start'} ],
					    abs($row[ $header{'chrom_stop'}] - 
						$row[ $header{'chrom_start'}] )];
    }
}

#end of docs
    
my $debug =0;

my $Datafile = '/projects/coccidiodies/variation/SnpReport_Cspp-3.txt';

my $Strategy = 'l'; # 'l' means locus based
                    # 'w' means window based
my $Window_size      = 1000; # only used if strategy='w'
my $Group_by_species = 1;
my $Min_seg_sites    = 5;
GetOptions(
	   'h|help' => sub { exec('perldoc', $0); 
			     exit(0);
			 },
	   'v|verbose!'   => \$debug,
	   'd|data:s'     => \$Datafile,
	   's|strategy:s' => \$Strategy,
	   'w|window:i'   => \$Window_size, # only used if strategy='w'
	   'g|groupbyspecies!' => \$Group_by_species,
	   );

$Datafile = shift @ARGV if @ARGV;

# validate the argument - our code is only 2 possible letters currently
$Strategy = lc(substr($Strategy,0,1));
if( $Strategy ne 'w' &&
    $Strategy ne 'l' ) {
    die("Can only understand 'w' (window) or 'l' (locus) strategy for this script\n");
}  

open( my $fh => $Datafile ) || die "cannot open $Datafile: $!";
my $i = 0;
chomp(my $hdr = <$fh>);
my %header = map { $_ => $i++ } split(/\t/,$hdr);

#debug: print the %header hash if you want to see what this is storing
#for my $h (  sort { $header{$a} <=> $header{$b} } keys %header ) {
#    print "$h $header{$h}\n";
#}


# A speedup to pre-index the by-species allele columns
my (@allele_cols,@sp_names, %sp_count, @HdrCols);

for my $col ( grep { /^a_(\S+)/ } keys %header ) {
    unless ( $col =~ /^a_(C[IP])/ ) {
	warn("unknown allele column '$col'\n");
    }
    $sp_count{$1}++;
    $sp_count{'ALL'}++; # for total number of sample count, when data is not partitioned by species
    push @allele_cols, [ $1, $header{$col} ];
}
@sp_names = sort grep { !/^ALL$/ } keys %sp_count;


if( $Strategy eq 'w' ) {
    push @HdrCols, qw(LINKAGE_GROUP CHROM_START WINDOW SIZE);
} elsif( $Strategy eq 'l' ) {
    push @HdrCols, qw(LINKAGE_GROUP CHROM_START LOCUS SIZE);
}

my @Methods = qw(tajima_D pi seg_sites);

if( $Group_by_species ) {
    for my $nm ( @sp_names ) {
	push @HdrCols, map { join(".",$_,$nm) } @Methods;
    }
} else {
    push @HdrCols, @Methods;    
}

print '#',join("\t", @HdrCols),"\n";
my %windows;

while(<$fh>) {
    chomp($_); # remove the trailing '\n'
    my @row = split(/\t/,$_); # split into an array

    # if we did re-sampling/removal of strains, we'd do it here
    my $allele_freq = {};
    for my $col ( @allele_cols ) {
	my $allele = $row[ $col->[1] ];
	next unless length($allele);
	$allele_freq->{$col->[0]}->{$allele}++;	
	$allele_freq->{'ALL'}->{$allele}++;  # for the all (Group_by_species = false) combined data	
    }
    
    my $window;
    if( $Strategy eq 'w' ) {
	 # process window-based accounting
	my $wbin = int($row[$header{'pos'}] / $Window_size);
	$window = sprintf("%s.%d",
			  $row[$header{'linkageGroup'}],
			  $wbin,
			  );

	unless( exists $windows{$window} ) {
	    $windows{$window} = { 'id'    => $wbin,
				  'pos'   => $wbin * $Window_size,
				  'length'=> $Window_size,
				  'sites' => []};
	}
	
    } elsif( $Strategy eq 'l' ) {
	 # process locus-based accounting
	my $g = $row[$header{'genes'}];	
	if( ! defined $g || length ($g) == 0 ) {
	    next; # skip this site because it doesn't fall within a defined locus
	}
	$g =~ s/\s+//g;
	$window = sprintf("%s.%s", $row[$header{'linkageGroup'}],
			  $g);
	unless( exists $windows{$window} ) {
	    $windows{$window} = { 'id'   => $g,
				  'pos'     => exists $genes{$g} ? $genes{$g}->[0] : '-1',
				  'length'  => exists $genes{$g} ? $genes{$g}->[1] : '-1',
				  'sites'=> []};
	}
    }
    push @{$windows{$window}->{sites}}, $allele_freq;
    last if $debug && (scalar keys %windows) > 5;
}
close($fh);

my $stats = Bio::PopGen::Statistics->new;
for my $window ( map { $_->[0] }
		 sort { $a->[1] <=> $b->[1] ||
			$a->[2] <=> $b->[2] } 
		 map { my $r = $_;		       
		       $r =~ s/[A-Z_]//g;
		       $r =~ s/,\S+//g;		       
		       [$_, split('\.',$r)] }
		 keys %windows ) {
    my ($lg,$wid) = split(/\./,$window);
    $lg = 'cimm_2.'.$lg;
    my ($id,$pos,$windowsize,$sites) = ( $windows{$window}->{'id'},
					 $windows{$window}->{'pos'},
					 $windows{$window}->{'length'},
					 $windows{$window}->{'sites'} );
    my $total_sites = scalar @{$sites};
    my @results;

    if( $Group_by_species ) {
	for my $sp ( @sp_names ) {
	    my ($seg_sites, $pi) = (0,0);	    
	    # calculate pi
	    for my $site ( @{$sites} ) {		
		my $ssh = 0;
		my @alleles = keys %{$site->{$sp}};
		my $sampsize = sum ( values %{$site->{$sp}});
		# check that alleles > 1 (not monomorphic in this population)
		# and that there are at more than 1 individual with a value here
		# (sampcount > 1)
		next unless( @alleles > 1 && $sampsize > 1);
		my $denom = $sampsize * ($sampsize - 1.0);	
		for my $a1 ( @alleles ) { 
		    $ssh += ( $site->{$sp}->{$a1} * ($site->{$sp}->{$a1} - 1)) / $denom;
		}
		$pi += 1.0 - $ssh;
		$seg_sites++;
	    }	    
	    my $D = 0;
	    if( $seg_sites >= $Min_seg_sites) {
		$D = $stats->tajima_D_counts($sp_count{$sp}, $seg_sites, $pi);
	    }
	    push @results, (sprintf("%.4f",$D),
			    sprintf("%.4f",$pi/$windowsize), 
			    $seg_sites);
	}
    } else {
	my $sp = 'ALL';
	my ($seg_sites, $pi) = (0,0);
	# calculate pi
	for my $site ( @{$sites} ) {	    
	    my $ssh = 0;
	    my @alleles = keys %{$site->{$sp}};
	    my $sampsize = sum ( values %{$site->{$sp}});
	    # check that alleles > 1 (not monomorphic in this population)
	    # and that there are at more than 1 individual with a value here
	    # (sampcount > 1)
	    next unless( @alleles > 1 && $sampsize > 1);
	    my $denom = $sampsize * ($sampsize - 1.0);	
	    for my $a1 ( @alleles ) { 
		$ssh += ( $site->{$sp}->{$a1} * ($site->{$sp}->{$a1} - 1)) / $denom;
	    }
	    $pi += 1.0 - $ssh;
	    $seg_sites++;
	}
	my $D = 0;
	if( $seg_sites >= $Min_seg_sites ) {
	    $D = $stats->tajima_D_counts($sp_count{$sp}, $seg_sites, $pi);
	}
	push @results, (sprintf("%.4f",$D),$pi/$windowsize, $seg_sites);
    }
    if( @results ) {
	print join("\t", $lg, $pos,$id,$windowsize, @results),"\n";
    }
}


