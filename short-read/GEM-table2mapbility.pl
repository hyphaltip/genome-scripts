#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use List::Util qw(sum);
my $genome = shift;
my $db = Bio::DB::Fasta->new($genome);
my $windowsize = 100;

my %hash;
my %data;
my $code;
my $seq;
my $readsize;
while(<>) {
    chomp;
    if( /^~([^~]+\S+)/) {
	$code = 'seq';
	$seq = $1;
    } elsif(/^~~(.+)/ ) {
	$code = $1;
    } elsif ( $code eq 'ENCODING' ) {
	my ($char,$rangel,$rangeh) = ( /\'(.)\'~\[(\d+)-(\d+)\]/ );
	$hash{$char} = [$rangel,$rangeh];
    } elsif( $code eq 'seq' ) {
	$data{$seq} .= $_;
    } elsif( $code eq 'READ LENGTH' ) {
	$readsize = $_;
    }
    
}

for my $seq ( keys %data) {
    open(my $fh => "| gzip -c > mapability/$seq.dat.gz") ||die $!;
    open(my $fhgc => "| gzip -c > gcWinds/$seq.dat.gz") ||die $!;
    my $len = length $data{$seq};
    for( my $i =0; $i < $len; $i+= $windowsize ){ 
	my $count = 0;
	my $j = 0;
	my @gc;
#	warn("i is $i\n");
	for( $j = 0; $j < $windowsize; $j++ ) {
	    last if $len <= ($j+$i);
	    my $ch = substr($data{$seq},($i+$j),1);
	    if( ! exists $hash{$ch} ) {
		die ("cannot find char '$ch' for $j ($i) len is $len\n");
	    }
	    my ($min,$max) = @{$hash{$ch}};
	    if( $min == 1 ) {
		$count++;
		my $read = $db->seq($seq,$i+$j => $i+$j+$readsize);
		if( $read ) {
		    push @gc, &gc_calc($read);
		}
	    }
	}
	printf $fh "%.2f\n",$count / $j;
	if( @gc ) {
	    printf $fhgc "%.12f\n", sum(@gc) / scalar @gc;
	} else {
	    printf $fhgc "NA\n";
	}
    }
    close($fh);
    close($fhgc);
}

sub gc_calc {
    my $str = shift @_;
    my ($gc) = ($str =~ tr/gcGC/gcGC/);
    return $gc / length($str);
}
