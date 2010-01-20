#!/usr/bin/perl -w
use strict;
use File::Spec;

my $min_length = 100_000;
my $ext = 'ss';
my $indir = File::Spec->rel2abs('ss');
my $keep  = File::Spec->rel2abs('ss/train');
my $skip  = File::Spec->rel2abs('ss/skip');

mkdir($keep);
mkdir($skip);
&cleanup($skip);
&cleanup($keep);

opendir(DIR,$indir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /\.\Q$ext\E$/;    
    my $infile = File::Spec->catfile($indir,$file);
    my $length = 0;
    open(my $fh => $infile ) || die $!;
    while(<$fh>) {
	if(/LENGTH\s+=\s+(\d+)/ ) {
	    $length = $1;
	    last;
	}
    }
    close($fh);

    if( $length >= $min_length ) {
	chdir($keep);	
	symlink(File::Spec->abs2rel($infile), 
		$file);
    } else {
	chdir($skip);
	symlink(File::Spec->abs2rel($infile), 
		$file);
    }
}


sub cleanup {
    my $dir = shift;
    opendir(DIR, $dir);
    for my $file ( readdir(DIR) ) {
	next unless $file =~ /\.\Q$ext\E$/;	
	unlink (File::Spec->catfile($dir,$file));
    }
}
