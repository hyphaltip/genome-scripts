#!/usr/bin/perl -w
use strict;

=head1 NAME 

fastatab_to_orthomcl  - turn FASTA tabular output into OrthoMCL BPO file

=head1 SYNOPSIS

 fastatab_to_orthomcl --score evalue|bit|bpr [-o outfile] inputfile1 inputfile2 ... 

=head1 DESCRIPTION

Comand line options:
  -s/--score    scaled|evalue|bit|bpr -- report scaled E-value,E-value, bitscore, or bits per residue
  -d/--dbsize             -- database size for scaling the E-value (default is 50,000)
  -o/--out                -- optional outputfile to write data, 
                             otherwise will write to STDOUT
  -h/--help               -- show this documentation

The expected input is tabular columns with the standard NCBI -m9 columns.

 queryname
 hit name
 percent identity
 alignment length
 number mismatches 
 number gaps
 query start  (if on rev-strand start > end)
 query end 
 hit start (if on rev-strand start > end)
 hit end 
 evalue
 fasta score

Plus the additional columns
 UNDEF
 bit score
 sw-score
 percent similar
 query length
 hit length

=head1 AUTHOR - Jason Stajich

Jason Stajich jason_at_bioperl.org

=cut


use strict;
use Getopt::Long;
my %data;

my ($score,$outfile,$debug) = qw(evalue);
my $dbsize = 50_000;
GetOptions(
	   'v|verbose!'          => \$debug,
	   's|score:s'           => \$score,
	   'o|out|outfile:s'     => \$outfile,
	   'd|db|dbsize:i'       => \$dbsize,
	   'h|help'              => sub { exec('perldoc',$0); exit; }
	   );

if( ! defined $outfile ) { 
    die("need and outfile with -o\n");
}
open(my $outfh, ">$outfile.bpo") || die("$outfile.bpo : $!");
open(my $ggfh,  ">$outfile.gg")  || die("$outfile.gg : $!");

my $counter = 1;
my %names;
while(<>) {
    my ($qname,$hname, $percent_id, $aln_len,
	$mm, $gaps, 
	$qstart, $qend,
	$hstart,$hend,
	$evalue, $fasta_score,
	$bitscore,undef, $sw_score,
	$percent_sim, 
	$qlength, $hlength, 
	$qgaps, $hgaps)  = split;
    my $score_value = $evalue;
    if( $score =~ /bpr/i ) {
	$score_value = $bitscore / abs($qend - $qstart);
    } elsif( $score =~ /bit/i ) {
	$score_value = $bitscore;
    } elsif( $score =~ /scaled|evaluescaled|escaled/i ) {
	$score_value = sprintf("%.2g",$dbsize * $qlength * 1 / 2**($bitscore));
    }
    for my $name ( $qname,$hname) { 
	my ($sp,$gname) = split(/\|/,$name);
	$names{$sp}->{$name} = 1;
    }
    # FROM ORTHOMCL DOCS
    # 1;At1g01190;535;At1g01190;535;0.0;97;1:1-535:1-535.
    # 2;At1g01190;535;At1g01280;510;2e-56;29;1:69-499:28-474.
    # 3;At1g01190;535;At1g11600;510;1e-45;27;1:59-531:21-509.
    
    # Each line represents each query-subject similarity
    # relation. And all the info is separated by ";", which
    # are, in order, similarity id, query id, query length,
    # subject id, subject length, BLAST E-value, percent
    # identity, HSP info (each HSP is in the format of
    # HSP_id:query_start-query_end:subject_start-subject_end.
    # different HSP info are seperated by "." )
    
    # IMPORTANT: 
    
    # BPO file is a parsing result from BLAST, so for each
    # query gene id, its hits can't be scattered in the file,
    # but should be listed in ajacent lines.  
    
    # only 1 HSP for SSEARCH/FASTA VALUES
    # so printing 1: for the HSP section

    my $hsp_info = sprintf("%d-%d:%d-%d",$qstart+1,$qend+1,$hstart,$hend);
    printf $outfh "%lu;%s;%d;%s;%d;%s;%d;1:%s.\n",
    ( $counter++,
      $qname,
      $qlength,
      $hname,
      $hlength,
      # defined $evalue ? $evalue : '',
      $score_value,
      $percent_id || 0,
      $hsp_info,
      );
    last if $debug && $counter > 1000;
}

for my $sp ( keys %names ) {
    printf $ggfh "%s: %s\n", $sp, join(" ",sort keys %{$names{$sp}});
}
