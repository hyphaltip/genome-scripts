#!/usr/bin/env perl
use strict;
use warnings;
use File::Path qw(make_path);
use Getopt::Long;
use Bio::SeqIO;

=head1 make_circos_karyotype

perl make_circos_karyotype.pl -i genome.fasta > 
=head1 USAGE


=cut

my ($in,$prefix,$out,$min_len);
my $max_colors = 360;
$min_len = 20_000;
GetOptions(
    'i|input:s' => \$in,
    'o|out:s'   => \$out,
    'm|min:s'   => \$min_len,
    'maxcolor:s' => \$max_colors,
    'p|prefix:s' => \$prefix);

$in ||= shift @ARGV;

my $seqio;
if( $in ) {
    $seqio = Bio::SeqIO->new(-format => 'fasta', -file => $in);
} else {
    $seqio = Bio::SeqIO->new(-format => 'fasta', -fh => \*STDIN);
}

my $ofh;
if( $out ) {
    open($ofh => ">$out") || die $!;
} else {
    $ofh = \*STDOUT;
}

my $i = 1;
while(my $seq = $seqio->next_seq ) {
    next if $seq->length < $min_len;
    my $id = $seq->display_id;
    my ($n);
    if( $id =~ /(\d+)$/) {
	$n = $1;
    } else { $n = $i; }
    my $color = sprintf("hue%03d",$i % $max_colors);
    if( $prefix ) {
	$id = $prefix ."_".$id;
    }
    printf $ofh "chr - %s scf%s %d %d %s\n", $id, $n, 0, $seq->length, $color;
    $i++;
}

