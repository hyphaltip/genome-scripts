#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;
use Bio::DB::SeqFeature::Store;
use Bio::Graphics::Wiggle::Loader;
use Env qw(GBROWSE_USER GBROWSE_PASSWORD);
die("must have GBROWSE_USER and GBROWSE_PASSWORD defined or fix this script\n") unless defined $GBROWSE_USER && defined $GBROWSE_PASSWORD;
my $mapfile;
my $wig;
my $wigdir = 'wigdb';
GetOptions(
	   'm|mapfile:s' => \$mapfile,	   
	   'w|wig|wiggle:s' => \$wig,
	   'wd|wigdb:s' => \$wigdir,
	   );

mkdir($wigdir) unless -d $wigdir;
my $wigdb = Bio::Graphics::Wiggle::Loader->new($wigdb);
open(my $fh => $wig) || die $!;
$wigdb->load($fh);
warn("$wigdb\n");

my $data = $wigdb->featurefile('featurefile');
warn("$data\n");

my @r = $data->values(1,10);
print "@r\n";
my $ref_genome_index = 0;
my @aln_lookups;

open(my $mapfh => $mapfile ) || die "$mapfile: $!";

while(<$mapfh>) {
    my ($aln_id, @line) = split;
    my ($chrom,$start,$end) = map { $line[$ref_genome_index + $_] } 0..2; 
    
    next if $chrom eq 'NA';
    $aln_lookups[$aln_id] = [$chrom,$start,$end];
}

my($user,$pass,$dbname,$host) = ($GBROWSE_USER,$GBROWSE_PASSWORD,
				 'gb_neurospora_crassa_cons',
				 'web.fungalgenomes.org');

my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );

my $segment = $dbh->segment('Ncra_OR74A_chrVI_contig7.4', 131800 => 134400);
my @scores = $segment->features('phastcons:phastCons_transcript-Intron');
for my $s ( @scores ) {
    print ref($s),"\n";
    print "$s\n";
    print $s->score,"\n";
}
