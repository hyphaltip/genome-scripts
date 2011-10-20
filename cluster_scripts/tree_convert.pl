#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;
use Getopt::Long;
my $outformat = 'nhx';
my $informat = 'newick';

GetOptions(
	   'if|i:s' => \$informat,
	   'of|o:s' => \$outformat);

my $in = Bio::TreeIO->new(-format => $informat,
			  -file   => shift);
my $out = Bio::TreeIO->new(-format => $outformat);
while( my $t = $in->next_tree) {
    $out->write_tree($t);
}
