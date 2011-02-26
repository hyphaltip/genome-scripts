#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;
my ($db,$output);
GetOptions(
	   'd|db:s'     => \$db,
	   'o|output:s' => \$output,
	   );
# simple start-stop features to seq
my $dbh = Bio::DB::Fasta->new($db || shift);
my $out;
if( $output) {
#    $out = Bio::SeqIO->new(-format => 'fasta',
#			   -file => ">$output");
    open($out => ">$output") || die "cannot open $output: $!\n";
} else {
#    $out = Bio::SeqIO->new(-format => 'fasta');
    $out = \*STDOUT;
}
while(<>) {
    next if /^\#/;
    chomp;
    my @row = split(/\t/,$_);
    my ($start,$end) = sort { $a <=> $b } ($row[3],$row[4]);
    ($start,$end) = ($end,$start) if $row[6] eq '-' || $row[6] eq '-1';
    my $seq = $dbh->seq($row[0], $start => $end );
    my ($name);
    if($row[-1] =~ /(ID|Parent|Gene)=([^;]+)/) {
	$name = $2;
    } else {
	$name = sprintf("%s_%d-%d",$row[0],$start,$end);
    }
#    $out->write_seq(Bio::Seq->new(-id => $name,
#				  -seq => $seq,
#				  -desc => sprintf("%s:%d..%d",
#						   $row[0],$start,$end)));
    printf $out ">%s %s:%d..%d length=%d score=%s\n%s\n",$name,$row[0],$start,$end,length($seq),$row[5],$seq;
}
