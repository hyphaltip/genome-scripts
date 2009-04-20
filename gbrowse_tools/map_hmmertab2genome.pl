#!/usr/bin/perl -w
use strict;
use Env qw(HOME);

=head1 NAME

 map_hmmertable2genome.pl - map Pfam domains from protein space to genome coordinates

=head1 SYNOPSIS

 map_hmmertable2genome.pl -u username -p password -db database -h dbhost -t HMMERtable > HMMERTable.genome_coordinates.gff3

=head1 DESCRIPTION

Map Pfam domains (or any protein feature) from protein coordinates to genomic coordinates.  Expects output from hmmer2table.pl 
which is in bioperl/scripts/searchio/hmmer_to_table.PLS

=head1 AUTHOR

Jason Stajich, jason-at-bioperl.org

=cut

use Bio::Coordinate::GeneMapper;
use Bio::DB::SeqFeature::Store;
use Bio::Location::Simple;

use Getopt::Long;
my ($user,$pass,$dbname,$host) = qw(gbUser gBr0wse);

$host ='localhost';
my $prefix;
my $debug = 0;

my $src = 'HMMER_Pfam';
my %domains;
my ($tabfile);
GetOptions(
	   'v|verbose!'     => \$debug,
	   'u|user:s'       => \$user,
	   'p|pass:s'       => \$pass,
	   'host:s'         => \$host,
	   'db|dbname:s'    => \$dbname,
	   't|tab:s'        => \$tabfile,
	   );
$tabfile = shift @ARGV unless defined $tabfile;

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor   => 'DBI::mysql',
                                          -dsn       => $dsn,
                                          -user      => $user,
                                          -password  => $pass,
                                          );
my $fh;
if( $tabfile =~ /\.gz$/ ) {
 open($fh => "zcat $tabfile | ") || die $!;
} else {
 open($fh,$tabfile) || die $!;
}
my %seen;
while(<$fh>) {
    my ($gene_name, $qstart,$qend, $domain, $hstart,$hend, 
	$score,$evalue) = split;
    $seen{"$gene_name.$domain"}++;
    my ($gene) = $dbh->get_features_by_name($gene_name);
    if( ! defined $gene ) {
	warn("cannot find $gene_name\n");
	next;
    }
    my @exons;
#    warn("gene name is $gene_name\n");
    for my $mRNA ( $gene->get_SeqFeatures ) {
	for my $cds ( $mRNA->CDS ) {
#	    warn($cds->to_FTstring, "\n");
	    push @exons, $cds;
	}    
	my $genemap = Bio::Coordinate::GeneMapper->new(-in    => 'peptide',
						       -out   => 'chr',
						       #-cds   => $mRNA,
						       -exons => \@exons,
						       );
	my $newloc = $genemap->map(Bio::Location::Simple->new
				   (-start => $qstart,
				    -end   => $qend,
				    ));

	for my $exon ( $newloc->each_Location ) {
	    print join("\t", 
		       $gene->seq_id, $src, qw(translated_nucleotide_match),
		       $exon->start,
		       $exon->end,
		       $evalue,
			$exon->strand > 0 ? '+' : '-',
		       '.',
		       sprintf('ID=%s__%s.%d;Name=%s',
			       $domain, $gene_name,
			       $seen{"$gene_name.$domain"}, 
			       $domain)),"\n";
	}
    }
}

sub read_cnf {
    my ($user,$pass) = @_;
    if( -f "$HOME/.my.cnf") {
        open(IN,"$HOME/.my.cnf");
        while(<IN>) {
            if(/user(name)?\s*=\s*(\S+)/ ) {
                $user = $2;
            } elsif( /pass(word)\s*=\s*(\S+)/ ) {
                $pass = $2;
            }
        }
        close(IN);
    }
    return ($user,$pass);
}

