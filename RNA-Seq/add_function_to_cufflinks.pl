#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my $cutoff_evalue = 1e-3;
my ($cufflinks,$function_fasta_file,$hmmer_file);
$cufflinks = 'isoform_exp.diff';
my $filter;
GetOptions(
    'c|cuff:s' => \$cufflinks,
    'f|function:s' => \$function_fasta_file,
    'h|hmmer:s'    => \$hmmer_file,
    'e|evalue:s'   => \$cutoff_evalue,
    'filter!'      => \$filter,
    );

die "need a cufflinks file" unless defined $cufflinks && -f $cufflinks;

die "need a domain or a functional fasta file" unless (
    ( defined $function_fasta_file && -f $function_fasta_file ) || 
    ( defined $hmmer_file && -f $hmmer_file) );

my %funct;
my @header_extras;
if( $hmmer_file ) {
    open(my $fh => $hmmer_file) || die "cannot open $hmmer_file: $!";
    while(<$fh>) {
	next if /^\#/;
	my ($domain,$dom_acc, $tlen, $qname, $qacc, $qlen,
	    $evalue,$score, $bias, $n,$ntot, 
	    $c_evalue, $i_evalue, $dom_score, $dom_bias,@rest) = split(/\s+/,$_,23);
	my $desc = pop @rest;
	chomp($desc);
	next if $i_evalue > $cutoff_evalue;
	$funct{$qname}->{HMMER_DOMAINS}->{$domain} = $desc;
    }
    push @header_extras, 'HMMER_DOMAINS';
}
if( $function_fasta_file ) {
    open(my $fh =>"grep '^>' $function_fasta_file |") || die "cannot open $hmmer_file: $!";
    while(<$fh>) {
	if(/^>(\S+)\s+protein\s+Name:\"([^\"]+)\" AED/) {
	    $funct{$1}->{FUNCT_DESC} = $2;
	} elsif(/^>(\S+)\s+(.+)/) {
	    my ($gene,$desc) = ($1,$2);
	    chomp($desc);
	    $funct{$gene}->{FUNCT_DESC} = $desc;	    
	} else {
	    warn("$_\n");
	}
    }
    push @header_extras, 'FUNCT_DESC';
}

open(my $fh => $cufflinks) || die "$cufflinks: $!";

my $header = <$fh>;
my $i =0;
my @header = split(/\s+/,$header);
my %h = map { $_ => $i++ } @header;
if( ! exists $h{'test_id'} ) {
    die("unexpected header, cannot find gene_id\n");
}

print join("\t", @header, @header_extras), "\n";
while(<$fh>) {
    my @row = split;    
    my $gene = $row[ $h{'test_id'} ];
    if( $filter ) {
	next unless $row[ $h{'significant'} ] eq 'yes';
    }
    my ($domains,$function);
    my @extra_cols;
    for my $xtra ( @header_extras ) {
	if( exists $funct{$gene}->{$xtra} ) {	
	    if( ref($funct{$gene}->{$xtra}) =~ /HASH/ ) {
		# deal with the array reference
		push @extra_cols, join(";",map { sprintf("%s=%s",$_,$funct{$gene}->{$xtra}->{$_}) }
						     sort keys %{$funct{$gene}->{$xtra}});
	    } else {
		push @extra_cols, $funct{$gene}->{$xtra};
	    }
	}
    }
    print join("\t", @row,@extra_cols),"\n";
}



