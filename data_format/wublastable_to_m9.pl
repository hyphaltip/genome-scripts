#!/usr/bin/perl -w
use strict;
while(<>) {
 chomp($_);
my @line= split(/\t/,$_);
next if/^#/;
 print join("\t", 
 	$line[0], $line[1],$line[10], $line[6],$line[9], $line[12]+$line[14],  $line[17],$line[18], $line[20], $line[21], $line[2], $line[5]),"\n";
}
