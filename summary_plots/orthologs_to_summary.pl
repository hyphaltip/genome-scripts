#!/usr/bin/perl -w
use strict;

my $method = 'orthomcl';
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
    $seen{$row[ $header{'gene'}]} = 0;
}
close($fh);

my (%cc_only, %orthologs,%paralogs);
# an orthomcl file  with 1st col is ortholog id, the rest of cols are genes
while(<>) {
    my ($id,@genes) = split;
    my @lbic = map { s/\.t1//; s/lbic\|//; $_ } grep { /^lbic\|/ } @genes;
    my @ccin = map { my (undef,$n) = split('\|',$_);
		     $n =~ s/T\d+$//;
		 $n; } grep { /^ccin\|/ } @genes;
    if( @genes != @lbic + @ccin ) {
	die("did not extract all the genes for $id, @genes\n");
    }
    next unless @ccin;
    for my $c ( @ccin ) {
	$seen{$c}++;
    }
    
    if( scalar @ccin > 1 && @lbic) {
	for my $c ( @ccin ) {
	    push @{$paralogs{$c}}, grep { $_ ne $c } @ccin,@lbic;
	}
    } elsif( scalar @ccin == scalar @genes ) {
	for my $c ( @ccin ) {
	    push @{$cc_only{$c}}, grep { $_ ne $c } @ccin, @lbic;
	}
    } else {
	for my $c ( @ccin ) {
	    push @{$orthologs{$c}}, @lbic;
	}
    }
}

open($fh => ">plot/coprinus_orphans.tab");
print $fh join("\t", '#gene',
	       qw(chrom chrom_start chrom_stop)),"\n";

for my $g ( sort grep { $seen{$_} == 0 } keys %seen ) {
    print $fh join("\t",$g,@{$gene{$g}}),"\n";
}

#for my $g ( sort keys %cc_orphans ) {
#    print $fh join("\t", $g, @{$gene{$g}} ),"\n";
#}

close($fh);

open($fh => ">plot/coprinus_orthologs.tab");
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

open($fh => ">plot/coprinus_paralogs.tab");
print $fh join("\t", '#gene',
	       qw(chrom chrom_start chrom_stop chrom_strand targets method)),"\n";

for my $g ( sort keys %paralogs ) {
    print $fh 
	join("\t", ( map { $gene_ext{$g}->[$header{$_}] } 
		     qw(gene chrom chrom_start chrom_stop strand)),	     
	     join(",",@{$paralogs{$g}}),
	     $method),"\n";
}

close($fh);

open($fh => ">plot/coprinus_only.tab");

print $fh join("\t", '#gene',
	       qw(chrom chrom_start chrom_stop chrom_strand targets methods)),"\n";

for my $g ( sort keys %cc_only ) {
    print $fh 
	join("\t", ( map { $gene_ext{$g}->[$header{$_}] } 
		     qw(gene chrom chrom_start chrom_stop strand)),	     
	     join(",",@{$cc_only{$g}}),
	     $method),"\n";

}

close($fh);
