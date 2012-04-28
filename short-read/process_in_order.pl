#!/usr/bin/perl -w
use strict;

my ($fq1,$fq2,$base) = @ARGV;


open(my $fh1 => $fq1) || die "$fq1: $!";
open(my $fh2 => $fq2) || die "$fq2: $!";


open(my $ofh1 => ">$base\_1.fq") || die "$base\_1 : $!";
open(my $ofh2 => ">$base\_2.fq") || die "$base\_2: $!";
open(my $ofh3 => ">$base\_unpair.fq") || die "$base\_unpair: $!";

my $last_seen;
my $i = 0;

my ($name,$record) = get_record($fh1);   
my ($name2,$record2) = get_record($fh2);

while(! eof($fh1) || ! eof($fh2) ) {
    if( ! $name && $name2 ) {
	print $ofh3 $record2;
	next;
    } elsif( $name && ! $name2 ) {
	print $ofh3 $record;
    } elsif( ! $name && ! $name2 ) {
	exit;
    }
    
    my ($FC,$row,$col,$tile) = split(/:/,$name);
    my ($FC_2,$row_2,$col_2,$tile_2) = split(/:/,$name2);
    
    if( $FC_2 != $FC ) {
	die("Flowcells don't match, stopping\n");
    }

    if( $row == $row_2 &&
	$col == $col_2 &&
	$tile == $tile_2 ) {
	print $ofh1 $record;
	print $ofh2 $record2;
	($name,$record) = get_record($fh1);
	($name2,$record2) = get_record($fh2);   
    }
    

    while( $row_2 > $row && $col_2 > $col && $tile_2 > $tile) {	
	print $ofh3 $record;
	($name,$record) = get_record($fh1);
	($FC,$row,$col,$tile) = split(/:/,$name);
    }
}


sub get_record {
    my ($fh) = shift;
    if( eof($fh) ) { 
	return (undef,undef);
    }
    my $record = <$fh>;
    my $name;
    if( $record =~ /^@(\S+)\#\d+\/([12])/ ||
	$record =~ /^@(\S+\.\d+)\s+\S+:\d+\/([12])/  ) {
	$name = $1;	
    } elsif( $record =~ /^@(\S+):[YN]/) {
	$name = $1;	
    } else{
	warn("cannot parse name out for $record");
	next;
    }

    for ( 0..2) {
	$record .= <$fh>;
    }
    ($name,$record);	
}
