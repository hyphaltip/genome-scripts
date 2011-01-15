#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat');
		  
# assume short read fastq with 4 lines per record!

my $size = 1_000_000;
my $basename;
my $outdir;
GetOptions('s|size:i'  => \$size,
	   'o|outdir:s' => \$outdir,
	   'b|basename:s' => \$basename,
	   );
my @files = @ARGV;

for my $file ( @files ) {
    $file = File::Spec->rel2abs($file);
    my ($vol,$dir,$fname) = File::Spec->splitpath($file);
    my $odir = $outdir || $dir;

    my $fh;
    if( $fname =~ /(\S+)\.(gz|bz2)$/) {
	open($fh => "$uncompress{$2} $file |") ||die $!;
	$fname = $1;
    } else {
	open($fh => $file) || die $!;
    }
    my @name = split(/\./,$fname);
    my ($ext) = pop @name;
    my $f = join(".",@name);    
    my $i = 0;
    open(my $ofh => sprintf(">$odir/%s.p%02d.%s",$f,$i++,$ext)) || die $!;
    my $n = 0;
    while(<$fh>) {
	unless( /^@/ ) {
	    chomp;
	    die("out of register, got line:\n$_\n  but expected line to start with '\@'\n");
	}
	print $ofh $_;
	for ( 0..2) {
	    $_ = <$fh>;
	    print $ofh $_;
	}
	if( ++$n >= $size ) {
	    close($ofh);
	    open($ofh => sprintf(">$odir/%s.p%02d.%s",$f,$i++,$ext)) || die $!;
	    $n= 0;
	}
    }
    close($ofh);
}



