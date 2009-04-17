#!/usr/bin/perl -w
use strict;
my $evalue_cutoff = 1e-2;

my $domainq = 'kinase';
my $genefile = 'plot/coprinus_gene_summary.tab';
open(my $fh => $genefile) || die $!;
my %genes;
my $header = <$fh>;
$header =~ s/^\#//;
my $i =0;
my %header = map { $_ => $i++ } split(/\t/,$header);
my (%gene,%gene_ext,%seen);

while(<$fh>) {
    my @row = split;
    $row[ $header{'gene'} ] =~ s/\.(\d+)$//;
    
    $gene{$row[ $header{'gene'} ]} = [$row[ $header{'chrom'} ],
				      $row[ $header{'chrom_start'} ],
				      $row[ $header{'chrom_stop'} ]
				      ];
    $gene_ext{$row[ $header{'gene'} ]} = [@row];

    $gene{$row[ $header{'gene'} ]} =~ s/\.\d+$//;
    $gene{$row[ $header{'gene'} ]} = [$row[ $header{'chrom'} ],
				      $row[ $header{'chrom_start'} ],
				      $row[ $header{'chrom_stop'} ]
				      ];    
    $gene_ext{$row[ $header{'gene'} ]} = [@row];
    $seen{$row[ $header{'gene'}]} = 0;
}
close($fh);

mkdir("plot/domains") unless -d "plot/domains";
my %geneswdomain;
while(<>) {
    my ($gene,$start,$end,$domain,$hstart,$hend,$score,$evalue) = split;

    $gene =~ s/^([^:|]+)[:|]//;
    $gene =~ s/T\d+$//;
    
    $geneswdomain{$domain}->{$gene}++ if $evalue <= $evalue_cutoff;
}

for my $domain ( keys %geneswdomain) {
    open($fh => ">plot/domains/$domain.dat") || die $!;
    print $fh join("\t", '#gene',
		   qw(chrom chrom_start chrom_stop domain count)),"\n";
    for my $g ( sort keys %{$geneswdomain{$domain}} ) {
	print $fh join("\t", ( map { $gene_ext{$g}->[$header{$_}] } 
			       qw(gene chrom chrom_start chrom_stop)),	     
		       $domain, $geneswdomain{$domain}->{$g}), "\n";
    }
}
