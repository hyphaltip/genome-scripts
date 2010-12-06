#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Bio::Range;
use Env qw(HOME);
my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $gene_source = 'gene:NC10_CALLGENES_FINAL_2'; # default to Neurospora data
my $padding = 200; # bp to look on either side of a gene
my $output;

GetOptions(
	   'v|verbose!' => \$debug,
	   'u|user:s' => \$user,
	   'p|pass:s' => \$pass,
	   'host:s'   => \$host,
	   'db|dbname:s' => \$dbname,
	   
	   'g|gene:s'   => \$gene_source,
	   'padding:i'  => \$padding,
	   
	   'output:s'  => \$output,
	   );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );

if( $output && $output ne '-' ) { 
    open($output => ">$output" ) || die $!;
} else {
    $output = \*STDOUT;
}

my $iter = $dbh->get_seq_stream(-type => $gene_source);
my (undef,$from) = split(/:/,$gene_source);
print $output join("\t", qw(CHROM GENE_LEFT MRNA_LEFT START_LEFT STOP_LEFT STRAND_LEFT),
		   qw(GENE_RIGHT MRNA_RIGHT START_RIGHT STOP_RIGHT STRAND_RIGHT 
		      OVERLAP STRAND_MATCH SINGLE_MATCH OVERLAP_TYPE)), "\n";
my %seen;
while( my $gene = $iter->next_seq ) {
    my $name = $gene->name;    
    for my $mRNA ( $gene->get_SeqFeatures('mRNA') ) {
	my @overlaps = grep { $_->name ne $name } $mRNA->segment->features(-type => $gene_source);
	for my $overlap ( @overlaps ) {
	    for my $overlap_mRNA ( $overlap->get_SeqFeatures('mRNA') ) {		
		next if $overlap_mRNA->name eq $mRNA->name;
		next if $seen{$overlap_mRNA->name}->{$mRNA->name}++ ||
		    $seen{$mRNA->name}->{$overlap_mRNA->name}++;
		# obtain the exact feature types overlapping
		my $range = Bio::Range->new(-start => $overlap->segment->start,
					    -end   => $overlap->segment->stop);		
		my ($istart,$istop,$istrand) = $range->intersection($mRNA);
		my $i_range = $dbh->segment($gene->seq_id, $istart => $istop);
		my %types;
		for my $f ( $i_range->features(qw(CDS three_prime_utr five_prime_utr)) ) {
		    $types{$f->primary_tag}++;
		}
# print out this overlapping
		print $output join("\t",
				   $gene->seq_id,
				   $name, $mRNA->name,
				   $mRNA->start,$mRNA->stop,$mRNA->strand,
				   $overlap->name, $overlap_mRNA->name,
				   $overlap_mRNA->start,$overlap_mRNA->end,
				   $overlap_mRNA->strand,
				   join("..",$istart, $istop),
				   $mRNA->strand == $overlap_mRNA->strand ? 'same' : 'opp',
				   @overlaps == 1 ? 'yes' : 'no',
				   join(",", sort keys %types),
			       ), "\n";	   
	    }
	}
    }
    last if keys %seen && $debug;
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
