#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::PrimarySeq;
use Getopt::Long;

my ($query,$ref,$output,$infile);
my $debug = 0;
GetOptions(
	   'v|verbose!'   => \$debug,
	   'q|query:s'    => \$query,
	   'r|ref:s'      => \$ref,
	   'o|output:s'   => \$output,
	   'i|in|input:s' => \$infile,
	   );

$infile ||= shift @ARGV;
die " need an infile " unless $infile;
if( ! $output ) {
    my $name = $infile;
    $name =~ s/\.showcoords//;
    $output = $name . ".blocks"; 
}
open(my $fh => $infile) || die "$infile: $!";
my $line = <$fh>;
($ref,$query) = split(/\s+/,$line) unless $query && $ref;

$line = <$fh>;
if( $line !~ /^NUCMER/) {
    die("expected a NUCMER show-coord file not '$line'\n");
}
my $qdb = Bio::DB::Fasta->new($query);
my $rdb = Bio::DB::Fasta->new($ref);

warn("query=$query ref=$ref\n") if $debug;
mkdir($output) unless -d "$output";
my $block = 0;
while(<$fh>) {
    if( /^\[/ || /^\s+$/) {
	next;
    } else {
	my ($refstart,$refend,$qstart,$qend,@row) = split;
	my ($qname,$refname) = (pop @row,pop @row);
	warn("qname,refname = $qname,$refname\n") if $debug;
	mkdir("$output/$block") unless -d "$output/$block";
	my $out_query = Bio::SeqIO->new(-format => 'fasta',
					-file   => ">$output/$block/query.fa");
	$out_query->write_seq(Bio::PrimarySeq->new
			      (-id => sprintf("%s_%d-%d",$qname,$qstart,$qend),
			       -seq => $qdb->seq($qname,$qstart => $qend),
			       ));
	
	my $out_ref = Bio::SeqIO->new(-format => 'fasta',
				      -file   => ">$output/$block/ref.fa");
	$out_ref->write_seq(Bio::PrimarySeq->new
			      (-id => sprintf("%s_%d-%d",$refname,$refstart,$refend),
			       -seq => $rdb->seq($refname,$refstart => $refend),
			       ));
	print join(",", $block,$qname,$qstart,$qend,$refname,$refstart,$refend),"\n";
	$block++;
    }
    last if $debug;
}

