#!/usr/bin/perl -w

=head1 load_novoalign 

=head2 USAGE
cd /project/solexa/small_RNA_fungi/NC0001/novocraft

load_novoalign.pl --project NC0001 --code 2006AAXX_4 -i neurospora_crassa_OR74A_7.2006AAXX_4.all.novoalign.gz --description "BCGSC sequenced smallRNA pooled conditions (2007-09 sample)" --genome neurospora_crassa_OR74A_7.fasta --organism='Neurospora_crassa_OR74A.7' 

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Env qw(HOME);
use DBI;
use lib "$HOME/projects/taylor_projects/sReadDB/lib";
use sReadDB;

use constant MAX_LEN => 36;
use constant COMMIT_INTERVAL => 10000;

my $DBH;

END {
    $DBH->disconnect() if $DBH;
}


my ($user,$pass,$dbname,$host);
$dbname = 'sReadDB';
$host   = 'localhost';
my ($debug,$input,
    $project,$project_description,$code,
    $genome,$genome_gff,
    $organism,$read_file);
my ($params,$app) = qw(novoalign novoalign);
my $sequence_format = 'fasta';
GetOptions(
	   'v|verbose!'     => \$debug,
	   'u|user:s'       => \$user,
	   'p|pass:s'       => \$pass,
	   'host:s'         => \$host,
	   'db|dbname:s'    => \$dbname,
	   'i|input:s'      => \$input,
	   'params:s'       => \$params,
	   'app:s'          => \$app,
	   'project:s'      => \$project,
	   'description:s'  => \$project_description,
	   'code:s'         => \$code,
	   'g|genome:s'     => \$genome,
	   'gff|genomegff:s'=> \$genome_gff,
	   'organism:s'     => \$organism,
	   'sf|seqformat:s' => \$sequence_format,
#	   'readfile:s'     => \$read_file,
	   );


unless(  defined $dbname ) {
    die("no dbname provided\n");
}

unless( defined $project ) {
    # need a project name for this
    die("Must provide a project name with --project\n");
}

($user,$pass) = &sReadDB::read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s:mysql_local_infile=1:mysql_server_prepare=1',$dbname,$host);
$DBH = DBI->connect($dsn,$user,$pass,
		    {RaiseError => 1,
		     AutoCommit => 0});

# some primary keys that will be used
my ($project_id,$method_id,$read_id);

# load the project id 
$project_id = &sReadDB::get_project_id($DBH,$project,$project_description,$code);
die if( $project_id < 0 );

warn("project id is $project_id for $project,$project_description,$code\n");

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

if( $input =~ /\.gz$/ ) {
    open($fh => "zcat $input | ") || die $!;
} else {
    open($fh,$input) || die $!;
}

while(<$fh>) {
    # just process the header first
    if( /^\#\s+(novoalign\s+.+)/ ) {
	$params = $1;
    } elsif (/^\#/ || /^\s+$/ ) {
	next;
    } else { 
	last;
    }
}
$method_id = &sReadDB::get_method_id($DBH,$app,$params);

my %seen;
my $q_read = $DBH->prepare_cached(qq{SELECT read_id FROM read_info
				     WHERE project_id = ? AND
				     machid = ?});
my $i_read = $DBH->prepare_cached(qq{INSERT INTO read_info 
				  (project_id, machid, seq, qual, 
				   trimmed, trim_5, trim_3) 
				  VALUES (?,?,?,?,?,?,?)});

my $i_read_loc = $DBH->prepare_cached(qq{INSERT INTO read_location 
					(method_id,chromosome_id,
					 read_id,startmatch,endmatch,strand,
					 score,mismatches,quality,location) 
					 VALUES (?,?,?,?,?,?,?,?,?,LineFromText(?))
				     });
my ($machid, $id,$seqstr,$qualstr, $status, $alignscore, 
    $alignqual, $chrom, $offset, $strand, $mismatches,$x,$y);
my $i =0 ;
while(<$fh>) {
    next if /^\#/;
    chomp;
    ($id,undef,$seqstr,$qualstr, $status, $alignscore, 
     $alignqual, $chrom, $offset, $strand, $x,$y, 
     $mismatches) = split(/\t/,$_);
    if( ! defined $qualstr || ! defined $seqstr ) {
	warn("$_\n");
	next;
    }
    $id =~ s/^[>@]//;
    # at some point I may go back in and re-update these to have
    # full (untrimmed) seq read.
    my $trimmed_read_len = length($seqstr);
    $q_read->execute($project_id,$code."_R".$id);
    my ($read_id);
    if( $q_read->rows ) {
	($read_id) = @{$q_read->fetchrow_arrayref};
    }
#    $q_read->finish;

    unless( $read_id ) {
	$i_read->execute($project_id,$code."_R".$id,$seqstr,
			 $qualstr, $trimmed_read_len < length($qualstr),
			 0, $trimmed_read_len);
	$read_id = $DBH->{mysql_insertid};
	if( $@ || ! $read_id ) {
	    warn($DBH->errstr, "\n");
	    $i_read->finish;
	    $DBH->rollback;
	    last;
	}
#	$i_read->finish;
    }
    next unless defined $chrom;
    $chrom =~ s/^>//;
    $i_read_loc->execute($method_id, $chrominfo->{$chrom}->[0],
			 $read_id,
			 $offset, $offset + $trimmed_read_len,
			 $strand eq 'F' ? 1 : 0,
			 $alignscore,$mismatches, $alignqual,
			 sprintf("LINESTRING(0 %d,0 %d)",
				 $offset, $offset+$trimmed_read_len));
    # novocraft doesn't include the full prefix so we add the code in
    if( $@ ) {
	warn($DBH->errstr,"\n");
	$DBH->rollback;
	last;
    }
    $DBH->commit if $i++ && ($i % COMMIT_INTERVAL) == 0;
    last if $debug;
}
$i_read_loc->finish;
$i_read->finish;
$q_read->finish;
$DBH->commit;
close($fh);

