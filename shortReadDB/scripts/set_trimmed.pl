#!/usr/bin/perl -w

=head1 set_trimmed

=head2 USAGE

Requires a FASTQ input file

 cd /project/solexa/small_RNA_fungi/NC0001/3prime-MAP

 load_reads.pl --project NC0001 --code 2006AAXX_4 --description "BCGSC sequenced smallRNA pooled conditions (2007-09 sample)" --organism='Neurospora_crassa_OR74A.7' --input 200M6AAXX_4.stripped.fastq 

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
my ($debug,$input,$code,$project, $project_description, $organism);
GetOptions(
	   'v|verbose!'     => \$debug,
	   'u|user:s'       => \$user,
	   'p|pass:s'       => \$pass,
	   'host:s'         => \$host,
	   'db|dbname:s'    => \$dbname,
	   'i|input:s'      => \$input,
	   'project:s'      => \$project,
	   'description:s'  => \$project_description,
	   'code:s'         => \$code,
	   'organism:s'     => \$organism,
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

# load the project id 
my $project_id = &sReadDB::get_project_id($DBH,$project,$project_description,$code);

die if( $project_id < 0 );

my $fh;

if( $input =~ /\.gz$/ ) {
    open($fh => "zcat $input | ") || die $!;
} else {
    open($fh,$input) || die $!;
}


my ($id,$seq,$space,$qual);
my $i = 0;
my $u_read = $DBH->prepare_cached(qq{UPDATE read_info SET adaptor_trimmed = ?,trim_3=?
					 WHERE project_id = ? AND machid = ?});

$DBH->do("LOCK TABLES read_info WRITE");
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
	$u_read->execute($len < MAX_LEN,$len,$project_id,$id);
    }; 
    if( $@ ) {
	warn($@,"\n");
    }
    last if $debug && $i++ > 250;
    
}
$DBH->do("UNLOCK TABLES");
close($fh);
$DBH->do("ANALYZE TABLE read_info");

