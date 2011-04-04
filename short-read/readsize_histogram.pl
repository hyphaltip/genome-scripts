#!/usr/bin/perl -w
use strict;

for my $file ( @ARGV ) {
    open( my $in => $file) || die $!;
    
    my %sizes;
    open(my $ofh => ">$file.size.dat") || die $!;
    my ($id,$seq);
    while(<$in>) {
	if( /^>(\S+)/) {
	    my $x = $1;
	    if( $id ) {		
		my $count = 1;
		if( $id =~ /(\d+)-(\d+)/ ) {
		    $count = $2;
		}
		$sizes{length($seq)} += $count;
	    }
	    $id = $x;
	    $seq = undef;
	} else {
	    chomp;
	    $seq .= $_;
	}
    }
    if( $id ) {		
	my $count = 1;
	if( $id =~ /(\d+)-(\d+)/ ) {
	    $count = $2;
	}
	$sizes{length($seq)} += $count;
    }
    
    print $ofh join("\t", qw(LENGTH SIZE)),"\n";
    for my $len ( sort { $a <=> $b } keys %sizes ) {
	print $ofh join("\t", $len, $sizes{$len}),"\n";
    }
    close($ofh);
    open($ofh => ">$file.R") || die $!;
    print $ofh <<EOF
dat <- read.table("$file.size.dat",header=T,sep="\t");
    
EOF
    ;
}

