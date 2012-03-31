#!/usr/bin/perl

$filenameA = $ARGV[0];
$filenameB = $ARGV[1];
$filenameOut = $ARGV[2];

open $FILEA, "zcat $filenameA | " || die $!;
open $FILEB, "zcat $filenameB | " || die $!;

open $OUTFILE, "> $filenameOut" || die $!;

while(<$FILEA>) {
	print $OUTFILE $_;
	$_ = <$FILEA>;
	print $OUTFILE $_; 
	$_ = <$FILEA>;
	print $OUTFILE $_; 
	$_ = <$FILEA>;
	print $OUTFILE $_; 

	$_ = <$FILEB>;
	print $OUTFILE $_; 
	$_ = <$FILEB>;
	print $OUTFILE $_;
	$_ = <$FILEB>;
	print $OUTFILE $_;
	$_ = <$FILEB>;
	print $OUTFILE $_;
}
