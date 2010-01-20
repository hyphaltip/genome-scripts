#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;

my $mapfile;
my $dir = 'phastcons_output';
my $outdir=  'wig';
my $base = 'neurospora_crassa_OR74A_7';
my $debug = 0;
my $name = 'phastcons';
my $desc = 'PhastCons Conservation';
GetOptions(
	'm|map|mapfile:s' => \$mapfile,
	'd|i|dir|in:s'    => \$dir,
	'o|out|output:s'  => \$outdir,
	'v|verbose!'      => \$debug,
	'b|base|basename:s'=> \$base,
	'n|name:s'        => \$name,
	'desc:s'          => \$desc,
	);
mkdir($outdir) unless -d $outdir;

open(my $mapfh => $mapfile) || die "need a mapfile: $!";

# hardcoded to just take the position from the 1st genome which is
# Ncrassa in this instance

# parsing the mercator map file
my @aln_lookups;
while(<$mapfh>) {
    my ($aln_id, @line) = split;
    $aln_lookups[$aln_id]  = [ $line[0], $line[1], $line[2], $line[3] ];
}
close($mapfh);
opendir(WIG,$dir) || die "cannot open $dir: $!";
my %groups;
for my $wigfile( readdir(WIG) ) {
    next unless ( $wigfile =~ /^(\d+)\.([^\.]+)\.wig/);
    my ($wig_id,$group) = ($1,$2);
    next if $aln_lookups[$wig_id]->[0] eq 'NA';
    my $fname = File::Spec->catfile($dir,$wigfile);
    next if ( -z $fname || -d $fname );
    push @{$groups{$group}}, [$wig_id,$fname];
}

for my $group ( keys %groups ) {
    open(my $fh => ">$outdir/$base.$group.wig") || die $!;	
    printf $fh "track type=wiggle_0 name=\"%s\" description=\"%s\"\n",$name,$desc;
    for my $fset ( sort { $a->[0] <=> $b->[0] } @{$groups{$group}} ) {
	my ($wig_id,$fname) = @$fset;
	open(my $ifh => $fname) || die "$fname: $!\n";
	while(<$ifh> ) { 
	    if( /^fixedStep/ ) {
		$_ = sprintf("fixedStep chrom=%s start=%d step=1\n", 
			     $aln_lookups[$wig_id]->[0],
			     $aln_lookups[$wig_id]->[1]);
	    }
	    print $fh $_;
	}
	close($ifh);
	last if $debug;
    }
    close($fh);
}
