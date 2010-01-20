#!/usr/bin/perl -w
use strict;

my $dir = shift || '.';
opendir(DIR, $dir) || die $!;
my %vals;
for my $file ( readdir(DIR) ) {
    if( $file =~ /consEntropy\.C_(\d\.\d+)-EL_(\d+)\.out/ ) {
	my ($coverage,$exp_len) = ($1,$2);
	my $id = join(",", $coverage,$exp_len);
	open(my $fh => "$dir/$file") || die $!;
	while(<$fh>) {
	    if( /PIT=L_min\*H=(-?\d+\.\d+)\s+bits/) {
		$vals{$id}->[0] = $1;		
		last;
	    }
	}
	close($fh);
	if( ! defined $vals{$id}->[0]) {
	    warn("cannot parse bits from $dir/$file\n");
	}

    } elsif( $file =~ /^C_(\d\.\d+)-EL_(\d+)$/ ) {
	my ($coverage,$exp_len) = ($1,$2);
	my $id = join(",", $coverage,$exp_len);	
	open(my $fh => "$dir/$file") || die $!;
	my $last;
	while(<$fh>) {
	    $last = $_;
	}
	if( $last =~ /cover\s+(\d+\.\d+)/ ) { 
	    $vals{$id}->[1] = $1;
	} else {
	    warn("cannot parse coverage from $dir/$file\n");
	}
	close($fh);
    }
}

# print in order sorted by X axis
print "#", join("\t",qw(BITS COVERAGE LABEL)),"\n";
for my $val ( sort { $vals{$a}->[0] <=> $vals{$b}->[0] } keys %vals){
    next unless @{$vals{$val}} == 2;
    print join("\t", @{$vals{$val}}, $val),"\n";
}
