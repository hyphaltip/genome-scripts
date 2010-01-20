#!/usr/bin/perl -w

# author(s) Jason Stajich jason_AT_bioperl.org
#           Ian Holmes    ihh_AT_berkeley.edu

# this script transforms xrate output run with 
#  -pp (posterior probability trace)
# and (optionally) with a mercator input map file
# into a summary file in a simple tabular format or 
#  WIG format - http://genome.ucsc.edu/goldenPath/help/wiggle.html

# This script is generally applicable, but is written for case where 
# an alignment column is assigned a rate class (there are N) by xrate.
# Basically the phastcons doppelganger model.

# Output score can be several options specified by --score (or -s) option
#  (0) rate:  xrate rate based on class identifier and posterior probability
#  (1) 1- rate
#  (2) rate - 1/N
#  (3) (rate - 1/N) * N/(N-1)
#  (4) 1 - (rate - 1/N) * N/(N-1)  (conservation score)

# these have the following ranges:
# (0) rate: [1/N .. 1]
# (1) [(N-1)/N .. 0]
# (2) [0 .. (N-1)/N]
# (3) [0 .. N]
# (4) [N .. 0] 

use strict;
use File::Spec;
use Getopt::Long;
use List::Util qw(max min);
use Stockholm; # from DART

my %compress = ('.gz' => 'zcat',
		'.bz2'=> 'bzcat');
my @scorename = qw(rate 1-rate rate-class cons_score 1-cons_score);
my $mapfile; # this is the mercator map file
my $ext = 'xrate.pp'; 
my $desc = 'xRate conservation 3-class model';
my $name = 'xRate_3class';
my $odir = 'wig';
my $base = 'neurospora_crassa_OR74A_7';
my $score_type = 0; # one of 0-4
my $debug = 0;
my $prec = 6;
GetOptions('m|map:s'     => \$mapfile,
	   'n|name:s'    => \$name,
	   'd|desc:s'    => \$desc,
	   'o|odir:s'    => \$odir,
	   's|score:i'   => \$score_type,
	   'p|precison:i' => \$prec,
	   'v|verbose|debug!' => \$debug,
	   );

unless( $score_type =~ /^[0-4]$/ ) {
    die("unknown score type must be withing [0-4] not '$score_type'.\n");
}

die "need a mapfile" unless $mapfile;
mkdir($odir) unless -d $odir;

open(my $mapfh => $mapfile) || die "need a mapfile: $!";

# hardcoded to just take the position from the 1st genome which is
# Ncrassa in this instance

# parsing the mercator map file
my @aln_lookups;
while(<$mapfh>) {
    my ($aln_id, @line) = split;
    $aln_lookups[$aln_id]  = [ $line[0], $line[1], $line[2], $line[3] ];
}
$name .= "_$scorename[$score_type]";
my $dir = shift || die "need a dir";
opendir(FILE, $dir ) || die $!;
my $ofile =  File::Spec->catfile($odir,"$base.$name.wig");
open(my $outfh =>">$ofile") || die $!;
printf $outfh "track type=wiggle_0 name=\"%s\" description=\"%s\"\n",
    $name,$desc;
for my $file ( readdir(FILE) ) {
    next unless $file =~ /^(\d+)\.(\S+\.)?\Q$ext\E(\.(bz2|gz))?$/;
    my $comp = $3;
    my $aln_id = $1;
    my $infile = File::Spec->catfile($dir,$file);
    next if -z $infile;
    warn($file,"\n");
    if( ! exists $aln_lookups[$aln_id] ) {
	warn("unknown alignment block: $aln_id\n");
	next;
    }
    # skip things do not map  to reference
    next if( $aln_lookups[$aln_id]->[0] eq 'NA' );

    my $fh;
    if( $comp ) {
	open($fh=> "$compress{$comp} $infile |") || die "$compress{$comp} $infile: $!";
    } else {
	open($fh=> "<$infile") || die "$infile: $!";	
    }
    my $stock = Stockholm->from_filehandle($fh);

    my $cols = $stock->columns;
    my @gff = @{$stock->gf_GFF};
    my @rmean = map (0, 1..$cols); # rate per column
    my $maxState = 0;
    for my $gff (@gff) {
	my @f = split /\t/, $gff;
	if ($f[4] == $cols) {
	    if ($f[8] =~ /lgPost=([\-\d\.\+\-eE]+)/) {
		my $lgPost = $1;
		if ($f[2] =~ /^S(\d+)$/) {
		    my $state = $1;
		    $maxState = $state if $state > $maxState;
		    # this is the rate
		    $rmean[$f[3]] += $state * 2**$lgPost;
		}
	    }
	}
    }    
    printf $outfh "fixedStep chrom=%s start=%s step=1\n",
    $aln_lookups[$aln_id]->[0], $aln_lookups[$aln_id]->[1];
    my $max_rate = max(map { $_ / $maxState } @rmean);
    if( $max_rate == 0 ){
	warn("no rates? max rate is 0\n");
    }
    for my $pos (1..$#rmean) {
	# rate, further scaled by max_rate so it fits between 0 and 1
	my $rate = ( $rmean[$pos] / $maxState) / $max_rate;
	my $score = $rate;
	if( $score_type == 1 ) {
	    $score = 1 - $rate;
	} elsif( $score_type == 2 ) {
	    $score = $rate - 1 / $maxState;
	} elsif( $score_type == 3 ) {
	    $score = ($rate - 1 / $maxState) * ( $maxState/ ( $maxState - 1));
	} elsif( $score_type == 4 ) {
	    $score = 1 - ($rate - 1 / $maxState) * ( $maxState/ ( $maxState - 1));
	}
	printf $outfh "%.$prec"."f\n", $score;
    }
    last if $debug;
}

