#!/usr/bin/perl -w
use strict;

=head1 load_gff

=head2 USAGE

load_gff.pl --organism Neurospora_crassa_OR74A.7 --input neurospora_crassa_OR74A_7.NCBI_PASA2.gff3 [--genome neurospora_crassa_OR74A_7.fasta]

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Env qw(HOME USER);
use DBI;
use lib "$HOME/projects/taylor_projects/sReadDB/lib";
use sReadDB;

my ($dbtype,$user,$pass,$dbname,$host);
$user   = $USER;
$dbtype = 'Pg';
$dbname = 'sReadDB';
$host   = 'localhost';
my ($debug,$input,$organism,$genome,$genome_gff);
my $sequence_format = 'fasta';

GetOptions(
	   'v|verbose!'     => \$debug,
	   'u|user:s'       => \$user,
	   'p|pass:s'       => \$pass,
	   'host:s'         => \$host,
	   'db|dbname:s'    => \$dbname,
	   'dbtype:s'       => \$dbtype,
	   'i|input:s'      => \$input,
	   'organism:s'     => \$organism,
	   'g|genome:s'     => \$genome,
	   'gf|genomegff:s' => \$genome_gff,
	   'sf|seqformat:s' => \$sequence_format,
	   );


unless(  defined $dbname ) {
    die("no dbname provided\n");
}

unless( defined $organism ) {
    # need a project name for this
    die("Must provide an organism name with --organism\n");
}
my $DBH = &sReadDB::connect($user,$pass,$host,$dbname,$dbtype);


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

if( scalar keys %{$chrominfo} == 0 ) {
    die("Must have pre-loaded the chromosome data");
}

my $fh;
if( $input =~ /\.gz$/ ) {
    open($fh => "zcat $input | ") || die $!;
} else {
    open($fh,$input) || die $!;
}

my $i_sth = $DBH->prepare_cached(qq{INSERT INTO feature
					(chromosome_id,fid,fname,
					 fparent,ftype,fsource,fstart,
					 fstop,fstrand,location)
				  VALUES (?,?,?,?,?,?,?,?,?,GeomFromText(?))});
my $lastline;
while(<$fh>) {
    next if /^\#/;
    chomp;
    my ($seqid,$src,$type,$start,$end,$score,$strand,$phase,$group) = split(/\t/,$_);
    next if $type eq 'gene';
    if( ! exists $chrominfo->{$seqid} ) {
	warn("no chrominfo for $seqid ($_)\n$lastline\n");
	next;
    }
    $lastline = $_;
    my %group = map { my @r = split(/=/,$_);
		      uc($r[0]) => [ split(/,/,$r[1]) ]; } split( /;/,$group);
    
    next if ( $src =~ /chromosome|scaffold|assembly/i || 
	      $type =~ /chromosome|scaffold|assembly/i );
    my ($ID,$name,$parent);
    if( exists $group{'ID'} ) {
	($ID) = @{$group{'ID'}};
    } else {
	warn("No Id for $lastline\n");	
    }
    if( exists $group{'NAME'} ) {
	($name) = ( @{$group{'NAME'}} );
	$name = unescape($name);
    } 
    if( exists $group{'PARENT'} ) {
	($parent) = @{$group{'PARENT'}};
    }
    eval { 
	my $strand_bit = $strand eq '+' ? 1 : 0;
	
	$i_sth->execute($chrominfo->{$seqid}->[0],
			$group{'ID'}->[0],
			$name,
			$parent,
			lc $type,
			lc $src,
			$start,
			$end,
			$strand_bit,
			sprintf("LINESTRING(%d %d,%d %d)",
				$chrominfo->{$seqid}->[0],$start, 
				$chrominfo->{$seqid}->[0],$end));
    };
    if( $@ ) {
	warn("$@ -- for \n");
	$DBH->rollback;
	last;
    }
}
close($fh);
$i_sth->finish;
$DBH->commit;
$DBH->disconnect;

sub unescape {
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}
