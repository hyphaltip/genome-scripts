#!/usr/bin/perl -w

=head1 load_SOAP

=head2 USAGE
cd /project/solexa/small_RNA_fungi/NC0001/SOAP

load_novoalign.pl --project NC0001 --code 2006AAXX_4 \
--input NC0001.Map_miRNA_stripped_3.SOAP.out.bz2
 --Description "BCGSC sequenced smallRNA pooled conditions (2007-09 sample)" \
 --genome neurospora_crassa_OR74A_7.fasta \
 --organism='Neurospora_crassa_OR74A.7' \
 --app SOAP \
 --params 'soap -a ../3prime-MAP/200M6AAXX_4.stripped.fastq -d neurospora_crassa_OR74A_7.fasta -r 2 -w 15 -c 41 -p 4 -s 8 -o NC0001.Map_miRNA_stripped_2.SOAP.out' \

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Env qw(HOME);
use DBI;
use File::Temp qw(tempfile);
my $tempbase = '/var/tmp';
use lib "$HOME/projects/taylor_projects/sReadDB/lib";
use sReadDB;

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat');

use constant MAX_LEN => 36;
use constant COMMIT_INTERVAL => 10000;
use constant COUNTER_THRESH  => 500_000;
my $DBH;

END {
    $DBH->disconnect() if $DBH;
}


my ($user,$pass,$dbname,$host,$dbtype);
$dbtype = 'Pg';
$dbname = 'sReadDB';
$host   = 'localhost';
my ($debug,$input,
    $project,$project_description,$code,
    $genome,$genome_gff,
    $organism,$read_file,$project_id);
my ($params,$app) = qw(soap SOAP);
my $sequence_format = 'fasta';
GetOptions(
	   'v|verbose|debug!' => \$debug,
	   'u|user:s'       => \$user,
	   'p|pass:s'       => \$pass,
	   'host:s'         => \$host,
	   'db|dbname:s'    => \$dbname,
	   'dbtype:s'       => \$dbtype,
	   'i|input:s'      => \$input,
	   'params:s'       => \$params,
	   'app:s'          => \$app,
	   'project:s'      => \$project,
	   'project_id:i'   => \$project_id,
	   'description:s'  => \$project_description,
	   'code:s'         => \$code,
	   'g|genome:s'     => \$genome,
	   'gff|genomegff:s'=> \$genome_gff,
	   'organism:s'     => \$organism,
	   'sf|seqformat:s' => \$sequence_format,
	   't|tmp:s'        => \$tempbase,
#	   'readfile:s'     => \$read_file,
	   );


unless(  defined $dbname ) {
    die("no dbname provided\n");
}

unless( defined $project || $project_id) {
    # need a project name for this
    die("Must provide a project name with --project or a known project id with --project_id\n");
}


$DBH = &sReadDB::connect($user,$pass,$host,$dbname,$dbtype);

# some primary keys that will be used
my ($method_id,$read_id);

# load the project id 
if( ! defined $project_id ) {
    $project_id = &sReadDB::get_project_id($DBH,$project,$project_description,$code);
}
die if( $project_id <= 0 );

my $fh;

if( $genome ) {
    if( &sReadDB::load_chromosomes_seq($DBH,$organism,$genome,$sequence_format) < 0 ) {
	die("could not load chromosome info");
    }
} elsif( $genome_gff ){
    die "not ready for this yet\n";
    if( &sReadDB::load_chromosomes_gff($DBH,$organism,$genome_gff) < 0 ) {
	die("could not load chromosome gff info");
    }
} 

my $chrominfo = &sReadDB::get_chromosome_info($DBH,$organism);

if( $input =~ /\.(gz|bz2)$/ ) {
    open($fh => "$uncompress{$1} $input | ") || die $!;
} else {
    open($fh,$input) || die $!;
}

$method_id = &sReadDB::get_method_id($DBH,$app,$params);

my $q_read = $DBH->prepare_cached(qq{SELECT read_id FROM read_info
				     WHERE project_id = ? AND
				     machid = ?});

my $i =0;


my ($tfh,$tempfile) = tempfile(DIR    => $tempbase, 
			       UNLINK => 1);
my $error = 0;
while(<$fh>) {
    next if /^\#/;
    last if $i++ > 100 && $debug;

    chomp;
    my ($read,$seqstr, $qualstr, $no_hits,
     undef, $trimmed_len,
     $strand, $chrom, $start,
     $mm,@mms) = split(/\t/,$_);

    if( ! defined $qualstr || ! defined $seqstr ) {
	warn("$_\n");
	next;
    }
    my ($read_id);
    $q_read->execute($project_id,$read);
    if( $q_read->rows ) {
	$q_read->bind_columns(\$read_id);
	$q_read->fetch;
    } else {
	warn("unknown read $read, (project=$project_id) did you already run load_reads?\n");
	$q_read->finish;
	$error = 1;
	last;
    }
    $q_read->finish;
    
    unless( defined $chrom && exists $chrominfo->{$chrom}) {
	warn("skipping $chrom, does not exist\n") if $debug;
	next;
    }
    $strand = ($strand eq '+' ? 1 : 0);    
#    $i_read_loc->execute($method_id, $chrominfo->{$chrom}->[0],		    
#			 $read_id, $start, $start + $trimmed_len,		    
#			 $strand, $mm,sprintf("LINESTRING(0 %d,0 %d)",
#			 $start, $start+$trimmed_len));
    if( $sReadDB::DBTYPE eq 'mysql' ) {
    print $tfh join("\t", '',$method_id, $chrominfo->{$chrom}->[0],		    
		    $read_id, $start, $start + $trimmed_len,		    
		    $strand,
		    '',$mm,'',),"\n";
	} else {
	print $tfh join("\t", $method_id, $chrominfo->{$chrom}->[0],
                    $read_id, $start, $start + $trimmed_len,
                    $strand,$mm),"\n";

}
#		    sprintf("LineFromText(LINESTRING(0 %d,0 %d))",
#			    $start, $start+$trimmed_len)),"\n";
    warn($i," records processed\n") if ($i > COUNTER_THRESH &&
					($i % COUNTER_THRESH) == 0);
}
close($fh);
close($tfh);
chmod(oct('0755'),$tempfile);

if( ! $error && $sReadDB::DBTYPE eq 'pg' ) {
	$DBH->do("COPY read_location (method_id,chromosome_id,read_id,startmatch,endmatch,strand,mismatches) FROM '$tempfile' DELIMITER '\t'");
	$DBH->do("UPDATE read_location set location= GeomFromText('LINESTRING(' || chromosome_id || ' ' || startmatch || ','|| chromosome_id || ' ' || endmatch || ')') WHERE location IS NULL");
	$DBH->commit;
} else {
$DBH->do("LOAD DATA INFILE '$tempfile' INTO TABLE read_location FIELDS TERMINATED BY '\t'");
}

