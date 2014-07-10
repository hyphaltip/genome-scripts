#!env perl
use strict;
use warnings;
use Getopt::Long;

my $hmmer_version = 3;
my $focus = 'query';
my $cutoff = 0.01;
GetOptions('version:i' => \$hmmer_version,
	   'evalue:s'  => \$cutoff,
	   'target'    => sub { $focus = 'target' },
	   'query'     => sub { $focus = 'query' },
	   'focus:s'   => \$focus,
    );
open(my $fh => shift || die ) || die $!;

if( $focus ne 'target' && $focus ne 'query' ) {
    die "focus must be target or query";
}
while(<$fh>) {
    next if /^\#/;
    my ($gene_name,$qstart,$qend,$domain,$hstart,$hend,$score,$evalue);
    my ($domacc,$tlen,$qacc,$qlen, $fullevalue,$fullscore,$fullbias,
	$n,$ntotal,$cvalue,$ivalue,$dombias,$envfrom,$envto,$acc,$desc);

    if( $hmmer_version == 2 ) {
	($gene_name,$qstart,$qend,$domain,$hstart,$hend,$score,$evalue) = split;
    } else {
	chomp;
	next if /^\#/ || /^\s+/;
	($domain,$domacc,$tlen,$gene_name,$qacc,$qlen,
	 $fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
	 $score,$dombias,
	 $hstart,$hend, $qstart,$qend,$envfrom,$envto,$acc,$desc) = 
	     split(/\s+/,$_,23);
	$evalue = $ivalue;
    }
    next if $evalue > $cutoff;
    
    if( $focus eq 'target' ) {
	print join ("\t", $domain,'HMMER','Domain',
		    $hstart,$hend, $score,
		    '.','.',sprintf("Target=%s %d %d;Pfam_Accession=%s;Description=%s",
				    $gene_name,
				    $qstart,$qend,
				    $domacc,
				    $desc)),"\n";
    } else {
	print join ("\t", $gene_name,'HMMER','Domain',
		    $qstart,$qend, $score,
		    '.','.',
		    sprintf("Target=%s %d %d;Pfam_Accession=%s;Description=%s",
			    $domain,
			    $hstart,$hend,
			    $domacc,
			    $desc)),"\n";

    }
}
