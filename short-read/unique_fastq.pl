use strict;
use warnings;
my $last = '';
my $count = 0;
while(<>) {
    my $id = $_;
    my $seq = <>;
    my $com = <>;
    my $qual = <>;
    if( $last eq $id ) {
	$count++;
    } else {
	print $id,$seq,$com,$qual;
    }	
    $last  = $id;
}

warn("$count duplicates\n");
