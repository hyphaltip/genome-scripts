#!/usr/bin/perl -w

# this is expecting bedtools coverage -counts output

for my $file ( @ARGV ) {
    my $ofile = $file;
    if( $ofile =~ s/\.txt/.norm.txt/ ) {
    } else {
	$ofile .= '.norm';
    }
    open(my $in => $file ) || die $!;
    open(my $out => ">$ofile") || die $!;
    my %data;
    while(<$in>) {
	my @row = split;
	push @{$data{$row[0]}}, [$row[0],$row[1],$row[2], sprintf("%.6f",$row[3] / abs($row[2]- $row[1])) ];
    }
    for my $chr ( sort keys %data ) {
	for my $r ( sort { $a->[1] <=> $b->[1] } @{$data{$chr}} ) {
	    print $out join("\t", @$r), "\n";
	}
    }
}
