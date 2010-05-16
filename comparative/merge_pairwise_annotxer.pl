#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $single_copy = 1; 
# if you don't want only single copy do
# --no-single


GetOptions('single!' => \$single_copy);
my @order = @ARGV;
my @pref = map { /vs-(\S+)\.tab/ } @order;

my $i =0;
my %dat;
my $ref_header;
my $ref_cols = 6;
my %ref_data;
for my $file ( @order ) {
    open(my $fh => $file) || die $!;
    my $header = <$fh>;
    $ref_header = $header unless defined $ref_header;
    while(<$fh>) {
	next if /NO_SYNTENIC_ALIGNMENT/ ||
	    /NO_GENES_IN_INTERVAL/;
	my @row = split;
	my $gene = $row[0];
	my @refd = splice(@row,0,$ref_cols);
	unless( defined $ref_data{$gene} ) {
	    $ref_data{$gene} = \@refd;
	}
	my $single = pop @row;
	# we skip the multi-orthologs (where a gene is split in the predictions)
	next if $single_copy && $single eq 'no';
	$dat{$refd[0]}->{$pref[$i]} = \@row;
    }
    $i++;
}

my @header = split(/\s+/,$ref_header);
my @final_header = splice(@header,0,$ref_cols);
for my $sp ( @pref ) {
    push @final_header,
    map { $sp . "_". $_ } @header;
}
print join("\t", @final_header),"\n";
for my $gene ( sort keys %dat ) {
    # we skip unless this is a single-copy gene found in all
    # species
    next unless scalar keys %{$dat{$gene}} == scalar @pref;
    print join("\t", @{$ref_data{$gene}},
	       map { @{$dat{$gene}->{$_}} } @pref),"\n";
}
