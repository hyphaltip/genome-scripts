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

--hmmer or --version can specify hmmer2 or hmmer3 
  '2' supposes you ran hmmer2table.pl to get your table (HMMER2 processing pipeline)
  '3' supposes you ran HMMER3 and are providing the file from --domtblout

=head1 AUTHOR

Jason Stajich, jason-at-bioperl.org

=cut

use Bio::Coordinate::GeneMapper;
use Bio::DB::SeqFeature::Store;
use Bio::Location::Simple;

use Getopt::Long;
my ($user,$pass,$dbname,$host) = qw();

$host ='localhost';
my $prefix;
my $debug = 0;

my $src = 'HMMER_Pfam';
my %domains;
my ($tabfile);
my $hmmer_version = 2;
my $cutoff = 0.01;
GetOptions(
	   'v|verbose!'     => \$debug,
	   'hmmer|version:i' => \$hmmer_version,
	   'u|user:s'       => \$user,
	   'c|cutoff:s'     => \$cutoff,
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
    my ($gene_name,$qstart,$qend,$domain,$hstart,$hend,$score,$evalue);
    if( $hmmer_version == 2 ) {
	($gene_name,$qstart,$qend,$domain,$hstart,$hend,$score,$evalue) = split;
    } else {
	my ($domacc,$tlen,$qacc,$qlen, $fullevalue,$fullscore,$fullbias,
	    $n,$ntotal,$cvalue,$ivalue,$dombias,$envfrom,$envto,$acc,$desc);

	chomp;
	next if /^\#/ || /^\s+/;
	($domain,$domacc,$tlen,$gene_name,$qacc,$qlen,
	 $fullevalue,$fullscore,$fullbias,$n,$ntotal,$cvalue,$ivalue,
	 $score,$dombias,
	 $hstart,$hend, $qstart,$qend,$envfrom,$envto,$acc,$desc) = split(/\s+/,$_,23);
	 $evalue = $ivalue;	 
    }
    next if $evalue > $cutoff;
    $gene_name =~ s/-mRNA-\d+$//;
    $seen{"$gene_name.$domain"}++;
    
    my ($gene) = $dbh->get_features_by_name($gene_name);
    unless ( defined $gene ) {
	($gene) = $dbh->features(-type => ['gene'],
				 -attributes => {locus_tag => $gene_name});
	unless( defined $gene ) {
	    warn("cannot find $gene_name\n");
	    next;
	}
    }

    my @exons;
    warn("gene name is $gene_name for $gene\n") if $debug;
    
    for my $mRNA ( $gene->get_SeqFeatures ) {
	for my $cds ( $mRNA->CDS ) {
	    warn($cds->to_FTstring, "\n") if $debug;
	    push @exons, $cds;
	}    
	if( ! @exons ) {
	    warn(" cannot find exons for $mRNA\n");
	    next;
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
		       sprintf('ID=%s__%s.%d;Name=%s;overlapping_gene=%s',
			       $domain, $gene_name,
			       $seen{"$gene_name.$domain"}, 
			       $domain,
			       $gene_name)),"\n";
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
