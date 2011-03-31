#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use POSIX qw(strftime);

my $querystr = "Fungi[Organism] and (internal transcribed spacer OR ITS OR ITS1 OR ITS2)";
my $date = strftime ("%Y%m%d",localtime);
my $ofile = "fungal_ITS_$date.txt";

my $query = Bio::DB::Query::GenBank->new(-db => 'nucleotide',
					 -query => $querystr,
					 );
print $query->count,"\n";

my $query2 = Bio::DB::Query::GenBank->new(-db => 'nucleotide',
					  -query => $querystr,
					  -maxids => $query->count);

my $count = $query2->count;

open(my $fh => ">$ofile")|| die $!;
print $fh join("\n", $query2->ids), "\n";
close($fh);
