#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;

use Getopt::Long;

my ($query,$ref,$output,$infile);

GetOptions(
	   'q|query:s' => \$query,
	   'r|ref:s'   => \$ref,
	   'o|output:s'=> \$output,
	   'i|in|input:s' => \$infile,
	   );
$infile ||= shift @ARGV;
open(my $fh => $infile) || die "$infile: $!";
my $line = <$fh>;
($query,$ref) = split(/\s+/,$line) unless $query && $ref;


$line = <$fh>;
if( $line !~ /^NUCMER/) {
    die("expected a NUCMER show-coord file not '$line'\n");
}
my $qdb = Bio::DB::Fasta->new($query);
my $rdb = Bio::DB::Fasta->new($ref);

while(<$fh>) {
    if( /^\[/) {
	next;
    } else {
	my ($qtart,$qend,$hstart,$hend,@row) = split;
	my ($qname,$hname) = (pop @row,pop @row);
	
	print join(",", $qname,$qtart,$qend,$hname,$hstart,$hend),"\n";
    }
}
