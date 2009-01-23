#!/usr/bin/perl -w

=head1 load_reads

=head2 USAGE

Requires a FASTQ input file

 cd /project/solexa/small_RNA_fungi/NC0001/raw

 load_reads.pl --project NC0001 --code 2006AAXX_4 --description "BCGSC sequenced smallRNA pooled conditions (2007-09 sample)" --organism='Neurospora_crassa_OR74A.7' --input 200m6aaxx_4.fastq --trimmed ../3prime-MAP/200M6AAXX_4.stripped.fastq

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Env qw(HOME USER);
use DBI;
use File::Temp qw(tempfile);
my $tempbase = '/var/tmp';
use lib "$HOME/projects/taylor_projects/sReadDB/lib";
use sReadDB;

use constant MAX_LEN => 36;
use constant COMMIT_INTERVAL => 10000;

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat');

my $DBH;

END {
    $DBH->disconnect() if $DBH;
}


my ($user,$pass,$dbname,$host,$dbtype);
$user   = $USER;
$dbtype = 'mysql';
$dbname = 'sReadDB';
$host   = 'localhost';
my ($debug,$input,$code,$project, $project_description, $organism,$trimmed);
GetOptions(
	   'v|verbose!'     => \$debug,
	   'u|user:s'       => \$user,
	   'p|pass:s'       => \$pass,
	   'host:s'         => \$host,
	   'db|dbname:s'    => \$dbname,
	   'dbtype:s'       => \$dbtype,
	   'i|input:s'      => \$input,
	   'n|name|project:s'      => \$project,
	   'description:s'  => \$project_description,
	   'code:s'         => \$code,
	   'organism:s'     => \$organism,
	   'trimmed:s'      => \$trimmed,
	   't|tmp:s'        => \$tempbase,
	   );


unless(  defined $dbname ) {
    die("no dbname provided\n");
}

unless( defined $project ) {
    # need a project name for this
    die("Must provide a project name with --project\n");
}

$DBH = &sReadDB::connect($user,$pass,$host,$dbname,$dbtype);

# load the project id 
my $project_id = &sReadDB::get_project_id($DBH,$project,$project_description,$code);

die if( $project_id < 0 );

my $fh;

if( $input =~ /\.(gz|bz2)$/ ) {
    open($fh => "$uncompress{$1} $input | ") || die $!;
} else {
    open($fh,$input) || die $!;
}


my ($tfh,$tempfile) = tempfile(DIR=> $tempbase, 
			       UNLINK => 1);
warn("tempfile is $tempfile\n");
my ($id,$seq,$space,$qual);
my $i = 0;
while(<$fh>) {
     chomp($id = $_);
     chomp($seq = <$fh>);     
     chomp($space = <$fh>);
     chomp($qual  = <$fh>);
     if( $space ne '+') {
	 warn("out of register?, space was $space, id was $id\n");
     }
     substr($id,0,1,'');
     if( $sReadDB::DBTYPE eq 'pg' ) {	 
	 print $tfh join("\t", $project_id, $id,$seq,$qual),"\n";
     } else {
	 # an empty default column
	 print $tfh join("\t", '',$project_id, $id,$seq,$qual),"\n";

     }
     last if $i++ > 100 && $debug;
     warn("$i records\n") if( $i > 500_000 && ($i % 500_000) ==0);
}
close($fh);
close($tfh);
chmod(oct('0755'),$tempfile);
if( $sReadDB::DBTYPE eq 'pg' ) {
     
    $DBH->do("COPY read_info (project_id, machid, seq, qual) FROM '$tempfile' DELIMITER '\t'"); 
} elsif( $sReadDB::DBTYPE eq 'mysql' ) {
    $DBH->do("LOAD DATA INFILE '$tempfile' INTO TABLE read_info FIELDS TERMINATED BY '\t'");
} else {
    warn("unknown dbtype $dbtype ($sReadDB::DBTYPE)\n");
}
if( $trimmed ) {
    if( $trimmed =~ /\.(gz|bz2)$/ ) {
	open($fh => "$uncompress{$1} $trimmed | ") || die $!;
    } else {
	open($fh,$trimmed) || die $!;
    }
    
    my $u_read = $DBH->prepare_cached(qq{UPDATE read_info SET adaptor_trimmed = ?,trim_3=?
					     WHERE project_id = ? AND machid = ?});
    
    $DBH->do("LOCK TABLES read_info WRITE") if $sReadDB::DBTYPE eq 'mysql';
    $i=0;
    while(<$fh>) {
	chomp($id = $_);
	chomp($seq = <$fh>);     
	chomp($space = <$fh>);
	chomp($qual  = <$fh>);
	if( $space ne '+') {
	    warn("out of register?, space was $space, id was $id\n");
	}
	substr($id,0,1,'');
	my $len = length($seq);
	eval {
	    $u_read->execute($len < MAX_LEN ? 1 : 0,$len,$project_id,$id);
	}; 
	if( $@ ) {
	    warn($@,"\n");
	}
	last if $debug && $i++ > 250;
    }
    $u_read->finish;
    if(  $sReadDB::DBTYPE eq 'mysql' ) {
	$DBH->do("UNLOCK TABLES");
	$DBH->do("ANALYZE TABLE read_info");
    } else {
	$DBH->commit;
    }
    close($fh);    
} else {
	$DBH->commit if ($sReadDB::DBTYPE eq 'pg');
}
