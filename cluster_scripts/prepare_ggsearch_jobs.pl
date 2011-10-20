#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

my $dir = shift;
opendir(DIR, $dir) || die $!;

my @dirs = shuffle grep { /^\d+$/ } readdir(DIR);

my $size = 100;
my $ct = 0;
my $i =0 ;
open(my $ofh => ">jobs/run_$dir.$i.sh") || die $!;
print $ofh "#!/bin/bash\n", "cd ~/blocks\n";
for my $d ( @dirs ) {
    print $ofh "/usr/local/bin/ggsearch36 -H -b 1 $dir/$d/query.fa $dir/$d/ref.fa > $dir/$d/align.GGSEARCH\n"; 
    if( ++$ct % $size == 0 ) {
	$i++;
	close($ofh);
	open($ofh => ">jobs/run_$dir.$i.sh") || die $!;
	print $ofh "#!/bin/bash\n", "cd ~/blocks\n";
    }
}
