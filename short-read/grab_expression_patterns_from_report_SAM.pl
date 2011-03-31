#!/usr/bin/perl -w
use strict;
use List::Util qw(sum max min);
use Statistics::Descriptive;
# expect a result from summarize_fwindow_report_SAM

my $cutoff = 5;
my $debug = 0;
open(my $debugfh => ">gene_classify.out") || die $!;
my (%patterns,%pattern_string);
my ($data_field) = 'TOTAL_NORM';

my ($file,$patfile) = @ARGV;

if( ! defined $file ) { die("must provide data file as 1st arg\n") }
open(my $fh=> $file ) || die("cannot open data file: $file\n");

if( ! defined $patfile ) { die("must provide patterns file as 2nd arg\n") }
open(my $dat => $patfile) || die ("cannot open patterns file: $patfile\n");

my $h = <$dat>;
if( $h !~ /^LABEL/ ) {
    die("unexpected format for patterns: $h\n");
}
my (undef,@h) = split(/\s+/,$h);
while(<$dat>) {
    my ($label,@pattern) = split;
    for my $nm ( @h ) {		
	$patterns{$label}->{$nm} = shift @pattern;
    }
}
close($dat);
my (@header,@rows);
while(<$fh>) {
    next if /^\#/;
    chomp;
    if( ! @header ) {
	@header = split(/\t/,$_);
    } else {
	push @rows, [split(/\t/,$_)];
    }
}

my $n = 0;
my %header = map { $_ => $n++ } @header;
my %sets;
my ($feature,$type,$chrom,$start,$end,$strand, @data) = @header;
my ($note) = pop @data;
$n = 6;
for my $r ( @data ) {
    if( $r =~ /^([^_]+_[^_]+)_(\S+)/) {
	$sets{$1}->{$2} = $n++;
    } else {
    	die("cannot parse header - stuck with $r\n");
    }
}
my $ii = 0;
my @data_sets = sort keys %sets;
for my $lbl ( keys %patterns ) {
    $pattern_string{ join("", map { $patterns{$lbl}->{$_} } @data_sets)} = $lbl;
}
my %gene_pattern;
my %gene2pattern;
for my $row ( @rows ) {
    my ( $fname) = ($row->[0]);
    my @data = map { $row->[$sets{$_}->{$data_field}] } @data_sets;


    my ($largest) = max @data;
    my ($smallest) = min @data;
    next if $largest < $cutoff;
    $smallest = 0.1 if $smallest == 0;
    
    my %results = map { $_ => '0' } @data_sets;
    if( &log2($largest / $smallest) > 1 ) {
	for my $seti ( @data_sets ) {
	    my $dati = $row->[ $sets{$seti}->{$data_field} ];	
	    $results{$seti} = '-',next if $dati == 0;
	    if( &log2($dati/$smallest) > 1 ) {
		$results{$seti} = '+';
	    } elsif( &log2($dati/$largest) < -1 ) {
		$results{$seti} = '-'; 
	    }
	}
    }   

    my $pat = join("", join("", map { $results{$_} } @data_sets));
    my $patid = 'NONE';
    if( exists $pattern_string{$pat} ) {
	$patid = $pattern_string{$pat};
	push @{ $gene_pattern{ $patid }}, [$fname,@data,$row->[-1]];
	$gene2pattern{$fname} = $pat;
    }
    print $debugfh join("\t", $fname,
			$pat, $patid,
			map { sprintf("%s (%s;%s)",$results{$_},
				      $row->[$sets{$_}->{$data_field}],
				      $_) } 		   
			@data_sets),"\n";
    last if $debug && $ii++ > 10;
}
while( my ($pattern,$grp) = each %gene_pattern ) {
    open(my $ofh => ">$pattern.dat") || die $!;    
    print $ofh join("\t", qw(FNAME), @data_sets, qw(NOTE)), "\n";
    print $ofh join("\n", map { join("\t", @$_) } @$grp),"\n";
}
sub log2 { log(shift) / log(2) }
