use strict;
use warnings;
my %dat;
while(<>) {
    chomp;
    my @row = split(/\t/,$_);
    my %group = map { split(/=/,$_) } split(';',pop @row);
    $row[2] = 'exon';
    my ($parent) = split(/\s+/,$group{Target});
    push @{$dat{$parent}}, 
    [@row, \%group];
}
for my $tran ( sort { my (undef,$at) = split(/_/,$a);
		      my (undef,$bt) = split(/_/,$b);
		      $at <=> $bt }
	       keys %dat ) {
    my $i = 1;
    	       
    my ($chr,$src,$min,$max,$strand);
    for my $exon ( @{$dat{$tran}} ) {
	$min = $exon->[3] if ! defined $min || $min > $exon->[3];
	$max = $exon->[4] if ! defined $max || $max < $exon->[4];
	$strand = $exon->[6] unless defined $strand;
	$src = $exon->[1] unless defined $src;
	$chr = $exon->[0] unless defined $chr;
    }
    print join("\t", $chr,$src,'gene',$min,$max,'.',$strand,'.',
	       sprintf("ID=%s.gene",$tran)),"\n";

    print join("\t", $chr,$src,'mRNA',$min,$max,'.',$strand,'.',
	       sprintf("ID=%s;Parent=%s.gene",$tran,$tran)),"\n";
    
    $strand = $strand eq '-' ? -1 : 1;
    for my $exon ( sort { $a->[3] * $strand <=> $b->[3] * $strand } 
		   @{$dat{$tran}} ) {
	my $group = pop @$exon;
	print join("\t",@$exon,
		   sprintf("ID=%s.e%d;Parent=%s;Target=%s",
			   $group->{ID},$i++,
			   $tran,
			   $group->{Target})), "\n";
    }
}

