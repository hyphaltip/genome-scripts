#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;
use Getopt::Long;

my ($format,$file,$out) = qw(fasta);
my ($matchgroup,$match,$source) = qw(match nucleotide_match);
my $percent_id = 60;
my $evalue_cutoff = 1e-10;
my $min_length = 200;
GetOptions(
	   'f|format:s'     => \$format,
	   'i|in|input:s'   => \$file,
	   'o|out|output:s' => \$out,
	   'mg|matchgroup:s'=> \$matchgroup,
	   'm|match:s'      => \$match,
	   's|source:s'     => \$source,
	   'e|evalue:s'     => \$evalue_cutoff,
	   'p|percentid:i'  => \$percent_id,
	   'ml|minlen:i'    => \$min_length,
	   );

$file = shift @ARGV unless defined $file;

if( $out ) {
    open($out => '>'.$out) || die $!;
} else{
    $out = \*STDOUT;
}

my $in = Bio::SearchIO->new(-format => $format,
			    -file   => $file);

print $out "#gff-version 3\n";
my ($match_id,$hspn) = (0,0);
while(my $r = $in->next_result ) {
    my $qname = $r->query_name;
    $source ||= $r->algorithm;    
    while( my $h = $r->next_hit ) {
	my $hname = $h->name;
	my @match; 
	while( my $hsp = $h->next_hsp ) {
	    last if $hsp->evalue > $evalue_cutoff;
	    push @match,
	    join("\t",
		 $hname,
		 $source,
		 $match,
		 $hsp->hit->start,
		 $hsp->hit->end,
		 $hsp->score,
		 $hsp->hit->strand > 0 ? '+' : '-',
		 '.', # frame
		 sprintf("ID=Align:match%05d;Parent=Align:%05d;Target=%s %d %d %s",
			 $hspn++,
			 $match_id,
			 $qname,
			 $hsp->query->start,
			 $hsp->query->end,
			 $hsp->query->strand < 0 ? '-' : '+')
		 );	    
	}
	if( @match ) {
	    print $out join("\t",
			    $hname, $source,$matchgroup,
			    $h->start('hit'),
			    $h->end('hit'),
			    $h->score,
			    $h->strand('hit') > 0 ? '+' : '-',
			    '.', # frame
			    sprintf("ID=Align:%05d;Name=%s;program=%s;programversion=%s;sourcename=%s;Target=%s %d %d %s",
				    $match_id++,
				    join("-",$qname,$hname,$r->algorithm),
				    $r->algorithm,
				    $r->algorithm_version,
				    $r->database_name,
				    $qname,
				    $h->start('query'),
				    $h->end('query'),
				    $h->strand('query') > 0 ? '+' : '-')
			    ),"\n";
	    print $out join("\n", @match), "\n";
	}
    }
}
