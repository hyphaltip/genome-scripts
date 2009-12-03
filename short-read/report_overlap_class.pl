#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Env qw(HOME);
use File::Spec;
use List::Util qw(sum);

use Bio::DB::SeqFeature::Store;

my %compress = ('bz2' => 'bzcat',
		'gz'  => 'zcat',
#		'dat'    => 'cat',
		'Z'   => 'zcat');
my ($user,$pass,$dbname); 
my $dir = 'feature_summary';
my $host = 'localhost';
my $debug = 0;
GetOptions(
	   'd|dir:s' => \$dir,
	   'db:s'    => \$dbname,
	   'u|user:s'=> \$user,
	   'p|pass:s'=> \$pass,
	   'host:s'  => \$host,
	   'v|verbose!'=> \$debug,
	   
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
mkdir($dir) unless -d $dir;

for my $file ( @ARGV ) {
    my (%ofh,%counts);
    my (undef,undef,$fname) = File::Spec->splitpath($file);
    my ($base,$ext) = split(/\.([^\.]+)$/,$fname,2);
    $base =~ s/\.dat$//;
    my $fh;
    if( exists $compress{$ext} ) {
	open($fh => "$compress{$ext} $file |") || die $!;
    } else {
	open($fh => "<$file") || die $!;
    }
    my $i = 0;
    while(<$fh>) {
	chomp;
	my ($name,$seq,$qual,$hitcount,
	    $flag, $length, $direction,
	    $chrom, $location_start, # 1 based
	    $types) = split(/\t/,$_);
	my $segment = $dbh->segment($chrom, $location_start => 
				    $location_start+$length);
	$direction = ($direction eq '+') ? 1 : -1;
	my $matched = 0;
	for my $f ( $segment->features() ) {
	    next if $f->type eq 'scaffold:chromosome'
		|| $f->type =~ /gene:/;
	    my $code = ($f->strand == $direction ) ? 'same' : 'diff';
	    $counts{$f->type}->{$code}++;	    
	    $matched ++;
	}
	unless( $matched ) { 
	    $counts{'nofeature'}->{same}++;
	}
	last if( $debug && $i > 1000);
	if( $i++ % 10_000 == 0 && $i > 1) {
	    warn("$base $i processed\n");
	}	
    }

    open(my $rpt => ">$dir/$base.classes_overlap.dat" ) || die $!;
    print $rpt join("\t", qw(CLASS COUNTS_SAME COUNTS_DIFF COUNTS_TOTAL PERCENT_SAME PERCENT_DIFF)),"\n";
    # print these in order of largest total number of reads that matches the
    # class using a Schwartizian transformation
    # sum( values %{$counts{$_}}) get the total number of 
    #  reads in each class (same or diff direction overlapping)
    for my $class ( map { $_->[0] }
		    sort { $b->[1] <=> $a->[1] }
		    map { [ $_, sum( values %{$counts{$_}})] }
		    keys %counts ) {
	my $sum = sum(values %{$counts{$class}});
	print $rpt join("\t", $class, 
			(map { $counts{$class}->{$_} || 0} qw(same diff)),
			$sum,
			(map { sprintf("%.1f%%",
				       (100 * $counts{$class}->{$_} || 0)/$sum) }
			 qw(same diff))),"\n";
    }
    close($rpt);
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
