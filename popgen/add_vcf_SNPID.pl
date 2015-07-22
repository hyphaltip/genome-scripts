use strict;
use warnings;
use Getopt::Long;

my $afumhack =0;

GetOptions('afumhack!' => \$afumhack);
while(<>) {
    unless( /^\#/ ) {
	my @row = split(/\t/,$_);
	my $chr = $row[0];
	if( $afumhack && $chr =~ /scf_(\d+)_(\S+)/ ) { # a hack for shortening some a.fumigatus contig/scf ids
	    $chr = 'scf'.int $1;
	}
	$row[2] = join(".",$chr,$row[1]);
	$_ = join("\t",@row);
    }
    print;
}

