#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $species = "Coprinus_cinereus";
GetOptions(
	's|species:s'=> \$species);
mkdir("Data_set");
mkdir("Data_set/$species");
mkdir("Data_set/$species/map");
mkdir("Data_set/$species/data");
my @lines; 

my ($gff,$gene2go) = @ARGV;
open(my $gff_fh => "<$gff" ) || die "$gff: $!";
open(my $genego_fh => "<$gene2go" ) || die "$gene2go: $!";

my %gene2chrom;
while(<$gff_fh>){ 
    next if /^\#/;
    my @line = split; 
    if ( $line[2] eq "mRNA" ) { 
	$line[-1] =~ s/ID=([^;]+);\S+/$1/; 
	push @lines, [$line[0], $line[3], $line[-1]];
	$gene2chrom{$line[-1]} = $line[0];
    } 
} 


my %gene2go;
while(<$genego_fh>){ 
    my ($gene,$go) = split;
    $gene2go{$gene}->{$go}++;
}


my $last;
my ($fh,$fh_data);
open(my $name => ">Data_set/$species/map_list") || die $!;
open(my $data => ">Data_set/$species/data_list") || die $!;

for my $l ( sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @lines ) { 
 if( ! defined $last || $l->[0] ne $last ) {
   close($fh) if defined $fh;
   close($fh_data) if defined $fh_data;
   print $name "Data_set/$species/map/".$l->[0].".map", "\n";
   open($fh, ">Data_set/$species/map/".$l->[0].".map") || die $!;
   print $data "Data_set/$species/data/".$l->[0].".data", "\n";
   open($fh_data, ">Data_set/$species/data/".$l->[0].".data") || die $!;
 }   
 print $fh join("\t", $l->[2], $l->[1], $l->[0]), "\n"; 
 if( exists $gene2go{$l->[2]} ) {
     for my $go ( keys %{$gene2go{$l->[2]}} ) {
	 print $fh_data join("\t", $l->[2], $go),"\n";
     }
 }
 $last = $l->[0];
}
close($name);
close($data);
close($fh);


