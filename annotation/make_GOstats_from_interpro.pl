#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;

=head1 NAME

make_GOstats_from_interpro - make a table for GOstats from an Interproscan TSV output

=head1 SYNOPSIS

perl make_GOstats_from_interpro.pl *.tsv > Organism.GO_table

=head1 DESCRIPTION

Generate a table suitable for processing with L<http://www.bioconductor.org/packages/release/bioc/html/GOstats.html> GOstats.
See L<http://www.bioconductor.org/packages/release/bioc/vignettes/GOstats/inst/doc/GOstatsForUnsupportedOrganisms.pdf> and 

=head1 AUTHOR 

 Jason Stajich E<lt>jason.stajich[at]ucr.eduE<gt>

=cut

my $GO_source = 'IEA';

GetOptions(
    'go|source:s' => \$GO_source,

    'h|help'      => sub { 
	exec('perldoc',$0);
	exit;
    }

    );
my %seen;

while(<>) {
    next if /^\#/;
    chomp;
    my ($gene,$crc64,$length,$analyis_method,
	$hit_dbid,$hitdesc,$qstart,$qend,$evalue,$hit_status,
	$date_run,$iprdomain,$iprdesc,
	$go_info) = split(/\t/,$_);

    if( $go_info ) {
	for my $go ( split(/\|/,$go_info) ) {
	    next if $seen{$go.$gene}++;
	    print join("\t", $go, $GO_source, $gene), "\n";   
	}
    }
}
