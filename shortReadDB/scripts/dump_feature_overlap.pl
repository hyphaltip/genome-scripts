#!/usr/bin/perl -w
use strict;

=head1 dump_feature_overlap

=head2 USAGE
 
 dump_featurer_overlap.pl [--project_id projectid] [--method_id methodid] [--organism Neurospora_crassa_OR74A.7] [--min 18] [--max 28] [--sep "\t"]

=head2 DESCRIPTION

--min     is the minimum read length to include in the report
--max     is the maximum read length to include in the report 

Values in [] indicate optional input. When real values are provided
they indicate the Defaults

=cut

use strict;
use Getopt::Long;
use Env qw(HOME);
use DBI;
use lib "$HOME/projects/taylor_projects/sReadDB/lib";
use sReadDB;

my $DBH;

END {
    $DBH->disconnect if defined $DBH;
}

my ($user,$pass,$dbname,$host,$prefix,$dbtype);
$dbname = 'sReadDB';
$host   = 'localhost';
$dbtype = 'Pg';
$prefix = '';
my ($debug,$min_size, $max_size,
    $organism,$only_trimmed,
    $project_id,$project, $code,
    $method_id, $app,$params);
#$min_size = 17;
#$max_size = 28;
my $sep = "\t";
$organism = 'Neurospora_crassa_OR74A.7';
my $dump_tRNA_only;
my $dump_rDNA_only;
my $no_overlaps;
GetOptions(
	   'v|verbose|debug!' => \$debug,
	   'u|user:s'         => \$user,
	   'pass|password:s'  => \$pass,
	   'host:s'           => \$host,
	   'db|dbname:s'      => \$dbname,
	   'dbtype:s'         => \$dbtype,

	   'p|o|output:s'     => \$prefix,	   
	   'min:i'          => \$min_size,
	   'max:i'          => \$max_size,
	   'sep:s'          => \$sep,
	   'tRNA!'          => \$dump_tRNA_only,
	   'rDNA!'          => \$dump_rDNA_only,
	   'n|none!'          => \$no_overlaps,
	   
	   'organism:s'     => \$organism,
	   'trimmed!'       => \$only_trimmed,

	   'method_id:i'    => \$method_id,
	   'app:s'          => \$app,
	   'params:s'       => \$params,

	   'project_id:i'   => \$project_id,
	   'project:s'      => \$project,
	   );


unless(  defined $dbname ) {
    die("no dbname provided\n");
}

unless( defined $organism || defined $organism ) {
    # need an organism name for this
    die("Must provide an organism name with --organism\n");
}

$DBH = &sReadDB::connect($user,$pass,$host,$dbname,$dbtype);
my $chrominfo = &sReadDB::get_chromosome_info($DBH,$organism);

my ($sth,@feature_types);
{ # capture feature types
    if( $dump_tRNA_only ) {
	$sth = $DBH->prepare(qq{
	    SELECT DISTINCT(ftype) FROM feature WHERE ftype = 'trna'});
    } elsif( $dump_rDNA_only) {
	$sth = $DBH->prepare(qq{
	    SELECT DISTINCT(ftype) FROM feature WHERE ftype = 'rdna'});
    } else {
	$sth = $DBH->prepare(qq{
	    SELECT DISTINCT(ftype) FROM feature WHERE ftype NOT IN ('rdna', 'trna')});
    }

    $sth->execute;
    my $ftype;
    $sth->bind_columns(\$ftype);
    while( $sth->fetch ){
 	push @feature_types, $ftype;
    }
    $sth->finish;
    warn( join("\n", reverse sort @feature_types),"\n") if $debug;
}
@feature_types = reverse sort @feature_types;
my ($read_id,$machid, $seq,$match_start, $match_end, $match_strand, 
    $trim_3,$is_trimmed,
    $feature_id, $fid, $fstart,$fend, $fstrand);

my ($method_str,$project_str,$trimmed) = ('','','');
if( $method_id ) {
    $method_str = "method_id = $method_id AND ";
}
if( $project_id ) {
    $project_str = "project_id = $project_id AND ";
}
if( $only_trimmed ) {
    if( $sReadDB::DBTYPE eq 'pg' ) {
	$trimmed = "ri.adaptor_trimmed = TRUE AND ";
    } else {
	$trimmed = "ri.adaptor_trimmed = 1 AND ";
    }
}
my $size_constrain_sql="";
if( $min_size && $max_size ) {
	$size_constrain_sql = "ri.trim_3 BETWEEN $min_size AND $max_size AND";
} elsif( $min_size ) {
	$size_constrain_sql = "ri.trim_3 >= $min_size AND";
} elsif( $max_size ) {
	$size_constrain_sql = "ri.trim_3 <= $max_size AND";
}
my ($sql);

if( $no_overlaps ) {
    if( $sReadDB::DBTYPE eq 'pg') {
	$sql = qq{SELECT ri.read_id, ri.seq, ri.machid,ri.trim_3,
		  ri.adaptor_trimmed,
		  r.startmatch,r.endmatch, r.strand
		  FROM read_info AS ri, read_location AS r
		      WHERE r.chromosome_id = ? AND $method_str $project_str 
                      $size_constrain_sql $trimmed 
		      ri.read_id = r.read_id AND
		      ri.read_id NOT IN  
		      (SELECT r.read_id FROM feature AS f 
		       WHERE r.chromosome_id = f.chromosome_id AND 
		       r.location && f.location)};
    } else {
	die("cannot do mysql right now for no-overlaps\n");
    }
    $sth = $DBH->prepare_cached($sql);
    open(my $ofh => ">$prefix\_NO_OVERLAPS.dat") || die $!;
    print $ofh join($sep, qw(READ_ID READ_SEQ READ_TRIM_LENGTH TRIMMED
			     CHROM MATCH_START MATCH_END 
			     MATCH_STRAND
			     )),
			     "\n";
    for my $cname (reverse sort { $chrominfo->{$a}->[0] <=> $chrominfo->{$b}->[0] }
		   keys  %{$chrominfo}  ) {
	warn("$cname\n") if $debug;
	my $wrote = 0;
	
	$sth->execute($chrominfo->{$cname}->[0]);
	$sth->bind_columns(\$read_id,\$seq,\$machid, \$trim_3, \$is_trimmed,
			   \$match_start, \$match_end, \$match_strand);
	while( $sth->fetch ) {
	    print $ofh join($sep, $machid, $seq,$trim_3, $is_trimmed,
			    $cname, $match_start, $match_end, 
			    $match_strand ? 1 : -1,
			    ),"\n";
	    $wrote++;
	}
	$sth->finish;
	$DBH->commit;
	last if $debug && $wrote;    
    }
    exit;
}

if( $sReadDB::DBTYPE eq 'pg' ) {
# I love you sub-queries    
    if( $dump_tRNA_only ) {
	$sql = qq{SELECT ri.read_id, ri.seq, ri.machid,ri.trim_3,ri.adaptor_trimmed,
		  r.startmatch,r.endmatch, r.strand,
		  f.feature_id, f.fid, f.fstart, f.fstop, f.fstrand
		      FROM read_info AS ri, read_location AS r, feature AS f
		      WHERE f.ftype = ? AND $method_str $project_str 
		      $size_constrain_sql $trimmed
		      r.chromosome_id = ? AND
		      r.chromosome_id = f.chromosome_id AND
		      ri.read_id = r.read_id AND r.location && f.location};
    } else {
	$sql = qq{SELECT ri.read_id, ri.seq, ri.machid,ri.trim_3,ri.adaptor_trimmed,
		  r.startmatch,r.endmatch, r.strand,
		  f.feature_id, f.fid, f.fstart, f.fstop, f.fstrand
		      FROM read_info AS ri, read_location AS r, feature AS f
		      WHERE f.ftype = ? AND $method_str $project_str 
		      $size_constrain_sql $trimmed
		      r.chromosome_id = ? AND
		      r.chromosome_id = f.chromosome_id AND
		      ri.read_id = r.read_id AND r.location && f.location
		      AND ri.read_id NOT IN (SELECT r.read_id FROM feature AS t 
					     WHERE t.ftype = 'trna' AND r.chromosome_id = t.chromosome_id AND 
					     r.location && t.location)};
    }

} else {
    if( $dump_tRNA_only ) {
	$sql = qq{SELECT ri.read_id, ri.seq, ri.machid,ri.trim_3,ri.adaptor_trimmed,
		  r.startmatch,r.endmatch, r.strand,
		  f.feature_id, f.fid, f.fstart, f.fstop, f.fstrand
		      FROM read_info AS ri, read_location AS r, feature AS f
		      WHERE f.ftype = ? AND $method_str $project_str 
		      $size_constrain_sql $trimmed
		      r.chromosome_id = ? AND
		      r.chromosome_id = f.chromosome_id AND
		      ri.read_id = r.read_id AND MBRIntersects(r.location,f.location)};
    } else {
	$sql = qq{SELECT ri.read_id, ri.seq, ri.machid,ri.trim_3,ri.adaptor_trimmed,
		  r.startmatch,r.endmatch, r.strand,
		  f.feature_id, f.fid, f.fstart, f.fstop, f.fstrand
		      FROM read_info AS ri, read_location AS r, feature AS f
		      WHERE f.ftype = ? AND $method_str $project_str 
		      $size_constrain_sql $trimmed
		      r.chromosome_id = ? AND
		      r.chromosome_id = f.chromosome_id AND
		      ri.read_id = r.read_id AND MBRIntersects(r.location,f.location) AND
		      ri.read_id NOT IN (SELECT r.read_id FROM feature AS t 
					 WHERE t.ftype = 'trna' AND r.chromosome_id = t.chromosome_id AND 
					 MBRIntersects(r.location, t.location))};
    }
    die "mysql not supported right now\n";
}

$sth = $DBH->prepare_cached($sql);
for my $ftype ( @feature_types ) {
    open(my $ofh => ">$prefix\_$ftype.dat") || die $!;
    print $ofh join($sep, qw(READ_ID READ_SEQ READ_TRIM_LENGTH TRIMMED
			     CHROM MATCH_START MATCH_END 
			     MATCH_STRAND
			     FEATURE FEATURE_START FEATURE_END FEATURE_STRAND)),
			     "\n";
    for my $cname (reverse sort { $chrominfo->{$a}->[0] <=> $chrominfo->{$b}->[0] }
		    keys  %{$chrominfo}  ) {
	warn("$cname\n") if $debug;
	my $wrote = 0;

	$sth->execute($ftype,$chrominfo->{$cname}->[0]);
	$sth->bind_columns(\$read_id,\$seq,\$machid, \$trim_3, \$is_trimmed,
			   \$match_start, \$match_end, \$match_strand,
			   \$feature_id, \$fid, \$fstart,\$fend, \$fstrand);
	while( $sth->fetch ) {
	    print $ofh join($sep, $machid, $seq,$trim_3, $is_trimmed,
			    $cname, $match_start, $match_end, $match_strand ? 1 : -1,
			    $fid, $fstart, $fend, $fstrand ? 1 : -1),"\n";
	    $wrote++;
	}
	$sth->finish;
	$DBH->commit;
	last if $debug && $wrote;    
    }
}
