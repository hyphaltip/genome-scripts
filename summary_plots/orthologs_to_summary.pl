#!/usr/bin/perl -w
use strict;

my $method = 'orthomcl';
my $genefile = 'plot/ci_gene_summary.tab';
open(my $fh => $genefile) || die $!;
my %genes;
my $header = <$fh>;
$header =~ s/^\#//;
my $i =0;
my %header = map { $_ => $i++ } split(/\t/,$header);
my (%gene,%gene_ext,%seen);

while(<$fh>) {
    my @row = split;
    $gene{$row[ $header{'gene'} ]} = [$row[ $header{'chrom'} ],
				      $row[ $header{'chrom_start'} ],
				      $row[ $header{'chrom_stop'} ]
				      ];
    $gene_ext{$row[ $header{'gene'} ]} = [@row];
    $seen{$row[ $header{'gene'}]} = 0;
}
close($fh);

my (%ci_orphans,%cocci_orphans, %orthologs);
# an orthomcl file  with 1st col is ortholog id, the rest of cols are genes
while(<>) {
    my ($id,@genes) = split;
#    my @cimm = map { s/\.t1/.2/; $_ } grep { /^CIMG_/ } @genes;
    my @cimm = map { s/\.t1//; $_ } grep { /^CIMG_/ } @genes;
    my @cpos = grep { /^cpos_/ } @genes;
    my @uree = grep { /^URET_/ } @genes;
    
    next unless @cimm;
    for my $c ( @cimm ) {
	$seen{$c}++;
    }
    # There are only Cimmitis genes
    if( scalar @cimm == scalar @genes ||
	(scalar @cpos + scalar @uree) == 0 ) {
	for my $c ( @cimm ) {
	    $ci_orphans{$c}++;
	}
    } elsif( (scalar @cpos + scalar @cimm) == scalar @genes ) {
	# There are only Cimmitis and Cposadasii genes
	for my $c ( @cimm ) {
	    push @{$cocci_orphans{$c}}, @cpos;
	}
    } else {
	for my $c ( @cimm ) {
	    push @{$orthologs{$c}}, @cpos,@uree;
	}
    }
}

open($fh => ">plot/ci_only.tab");
print $fh join("\t", '#gene',
	       qw(chrom chrom_start chrom_stop)),"\n";

for my $g ( sort keys %ci_orphans ) {
    print $fh join("\t", $g,@{$gene{$g}}),"\n";
}

close($fh);

open($fh => ">plot/ci_orphans.tab");
print $fh join("\t", '#gene',
	       qw(chrom chrom_start chrom_stop)),"\n";
for my $g ( sort grep { $seen{$_} == 0 } keys %seen ) {
    print $fh join("\t", $g, @{$gene{$g}}),"\n";
}
close($fh);

open($fh => ">plot/cocci_orphans.tab");
print $fh join("\t", '#gene',
	       qw(chrom chrom_start chrom_stop)),"\n";

for my $g ( sort keys %cocci_orphans ) {
    print $fh join("\t", $g,@{$gene{$g}}),"\n";
}

close($fh);

open($fh => ">plot/orthologs.tab");
print $fh join("\t", '#gene',
	       qw(src src_start src_stop src_strand 
		  target method)),"\n";
for my $g ( sort keys %orthologs ) {
    print $fh 
	join("\t", ( map { $gene_ext{$g}->[$header{$_}] } 
		     qw(gene chrom chrom_start chrom_stop strand)),
	     
		   join(",",@{$orthologs{$g}}),
	     $method),"\n";
}
close($fh);
