use strict;
use warnings;
print join("\t", qw(Strain Total Mapped Duplicate Unique)),"\n";
for my $file ( @ARGV ) {
    next unless ($file =~ /\.flagstat/);
    my ($base) = split(/\./,$file);
    open(my $fh => $file) || die $!;
    
    my ($total,$duplicate,$mapped);
    while(<$fh>) {
	if(/^(\d+).+ in total/ ) {
	    $total = $1;
	} elsif( /^(\d+).+ duplicates/) {
	    $duplicate = $1;
	} elsif( /^(\d+)\s+\+\s+\d+\s+mapped/) {
	    $mapped = $1;
	} else {
#	    warn($_);
	}
    }
    print join("\t", $base,$total,$mapped,$duplicate,$mapped - $duplicate),"\n";
	
}
