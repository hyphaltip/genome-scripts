#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $write_skip;
GetOptions(
	   'debug|skip!' => \$write_skip,
	   );

if( $write_skip  ) {
    open($skipfh => ">skip") || die $!;
}
while(my $id = <>) {
    my $seq = <>;
    my $desc =<>;
    my $qual = <>;
    
    if (($id =~ /\S+\s+\d+\:([YN]):\d+/) && $1 eq 'N') {
	print $skipfh $id,$seq,$desc,$qual if $write_skip;
	next;
    }
    
    print $id,$seq,$desc,$qual;
}
