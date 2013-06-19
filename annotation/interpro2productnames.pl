#!/usr/bin/perl -w
use strict;
use warnings;

#TODO - rework this to preferentiall use Panther families, but it will require a lookup to get the full name
# and fix the capitalization

=head1 NAME

interpro2productnames - generate product names for GenBank submission from an InterPro report

=head1 SYNOPSIS

interpro2productnames -i interpro.tab -o product_names.out

=head1 DESCRIPTION

Convert InterPro (tabular) results into product names for genes.

=head1 AUTHOR

Jason Stajich jason.stajich[at]ucr.edu

=cut

use Getopt::Long;

my ($interpro,$output);
my $debug = 0;
GetOptions(
    'i|ipr|interpro:s'   => \$interpro,
    'o|output:s'         => \$output,
    'v|verbose|debug:s'  => \$debug,
    );


my $fh;
if( $interpro ) {
    open($fh => $interpro) || die "cannot open $interpro: $!";
} else {
    $fh = \*ARGV;
}

my $ofh;
if( $output ) {
    open($ofh => ">$output") || die "cannot open $output for writing: $!";
} else {
    $ofh = \*STDOUT;
}

my %skip = map { $_ => 1} qw(IPR017870);

print $ofh join("\t", qw(PROTEIN PRODUCT DOMAIN_COUNT EC_NUM KEGG_ACC
METACYC_ACC REACTOME_ACC UNIPATHWAY)),"\n";
my %domains;
while(<$fh>) {
    chomp;    
    my ($protein,$crc64,$length,$analysis_method,
	$hit_dbid,$hitdesc,$qstart,$qend,$evalue,$hit_status,
	$date_run,$iprdomain,$iprdesc,
	$go_info,$pathways) = split(/\t/,$_);
    next if( defined $evalue && $evalue ne '-' && $evalue > 1e-20);
    push @{$domains{$protein}->{$analysis_method}}, [$hit_dbid, 
						     $hitdesc,
						     $qstart,$qend,$evalue,
						     $hit_status,
						     $iprdomain, $iprdesc,
						     $go_info,$pathways];
}

for my $protein (sort keys %domains)  {
    my %domain_set;
    my %pathways;
    #my $override_product;
    my $tm_count = scalar @{$domains{$protein}->{TMHMM} || []};
    for my $method ( keys %{$domains{$protein}} ) {
	next if $method =~ /Phobius|TMHMM/;
	next if ($method eq 'SUPERFAMILY' && exists $domains{$protein}->{Pfam});

	for my $domain ( @{$domains{$protein}->{$method}} ) {
	    my ($dbid,$desc,$qs,$qe,$evalue,$status,
		$iprdomain,$iprdesc,$go,$pathways) = @$domain;
	    my $len = ($qe - $qs);
	    next if defined $iprdomain && $skip{$iprdomain}; 

	    if( $pathways ) {
		for my $pthwy ( split(/\|/,$pathways ) ) { 
		    if( $pthwy =~ /(\S+):\s+(\S+)/ ) {
			my ($type,$val) = ($1,$2);
			if( $type eq 'KEGG' ) {
#			    $override_product = $iprdesc;
			    if( $val =~ /\+/ ) {				
				my ($keggid,$ec) = split(/\+/,$val);
				$pathways{$type}->{$keggid}++;
				$pathways{'EC'}->{$ec}++;
			    } else {
				warn("KEGG fmt ID is unexpected: $val");
				$pathways{$type}->{$val}++;
			    }
			} else { 
			    $pathways{$type}->{$val}++;
			}
		    }
		}
	    }
	    
	    if( ! $iprdesc && $method eq 'Pfam' ) {
		$domain_set{$desc}+= $len;
	    } elsif( $iprdesc ) {
		$domain_set{$iprdesc} += $len;
	    }
	}
    }
    my %seen;
    my $domain_count = scalar keys %domain_set;
    for my $dom ( sort { $domain_set{$b} <=> $domain_set{$a} }
		  keys %domain_set ) {
	my $val = $domain_set{$dom};
	next if $dom =~ /\S+[\s\-]+repeat(s)?(-like)?-containing/;
#	next if $dom =~ /\S+ repeats/;
	next if $dom =~ /helix-turn/;
	next if $dom =~ /repeat/;
	$dom =~ s/,\s+transmembrane//;
	next if $dom =~ /^ABC transporter(-like)?$/;
	$dom =~ s/\,\s*conserved site//;
	$dom =~ s/repeat signature//;
	$dom =~ s/\s*\,?\s*[CN]\-terminal\s*//g;
	$dom =~ s/,?\s*domain//;
	$dom =~ s/Domain of unknown function\s+//g;
	$dom =~ s/,\s*core(\s+domain)?//;
#	$dom =~ s/,\s*(conserved site|insertion protein)//g;
#	$dom =~ s/\s+cluster insertion//g;	
	next if $dom =~ /truncated/;
#	$dom =~ s/\s+homology/-like/;
	$dom =~ s/Polyadenylate binding protein, human types.+/Polyadenylate binding protein/;
	$dom =~ s/, active site//;
	$dom =~ s/,\s+(\S+)-fold//;
	$dom =~ s/_/-/g;
	$seen{$dom} += $val;
    }
    if( $debug ) {
	for my $dom ( sort { $seen{$b} <=> $seen{$a} } keys %seen ) {
	    warn("$dom $seen{$dom}\n");
	}
    }
    my @final = sort { $seen{$b} <=> $seen{$a} } keys %seen;
    if( @final > 1 ) {
	@final = ($final[1]); #,$final[2]);
    }
    
    $domain_count = scalar @final;
    my $product_name;
    if( ! $domain_count ) {
	$product_name = 'hypothetical protein';
#    } elsif( $override_product )  {
#	$product_name = $override_product;
    } else {
	$product_name = join(", ", @final);	
	$product_name .= sprintf(" domain-containing protein");
    }

    if( $product_name =~ /Nucleoporin, Nup133\/Nup155-like/) {
	$product_name = 'Nucleoporin, Nup133/Nup155-like';
    } 
    print $ofh join("\t", $protein, $product_name, $domain_count,
		    map { join(",", sort keys %{$pathways{$_} || {}}) }
		    qw(EC KEGG MetaCyc Reactome UniPathway)),"\n";
}
