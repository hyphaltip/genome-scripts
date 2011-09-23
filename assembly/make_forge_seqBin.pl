#!/usr/bin/perl -w

my ($reads_454_start, $reads_454_end,$reads_solexa_start,$reads_solexa_end) = @ARGV;


for(my $i = $reads_454_start; $i< $reads_454_end; $i++) {
 print join(" ", $i, '454','0'), "\n";
}

for(my $i = $reads_solexa_start; $i< $reads_solexa_end; $i++) {
 print join(" ", $i, 'solexa','0'), "\n";
}
