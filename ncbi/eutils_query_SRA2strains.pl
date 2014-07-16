#!env perl
use strict;
use warnings;

# SRA accession to Strain ID where it is encoded in a specific way (using the accessions from BioProject: PRJEB4236 in this case)
# provide as input a file with each SRA accession number on a line
# will print out the accession and then a strain designation if the metadata are encoded in the same way as this study
use Data::Dumper;
use File::Path;
use Getopt::Long;
use File::Spec;
use Text::CSV_XS qw(csv);
use LWP::Simple;
use XML::Simple;
use Encode;
use Cache::File;
use Bio::SeqIO;
use IO::String;

my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $SLEEP_TIME = 2;
my $cache_dir = "eutils_".$ENV{USER}.".cache";
my $cache_filehandle;
my $cache_keep_time = '2 day';

my $force = 0;
my $debug = 0;
my $retmax = 1000;
my $runonce = 0;
my $use_cache = 1;


my $DEBUG = 0;
GetOptions(
    'runonce!'  => \$runonce,
    'retmax:i'  => \$retmax,
    'f|force!'  => \$force, # force downloads even if file exists
    'cache!'    => \$use_cache,
    'v|d|debug'            => \$DEBUG,
    'f|force!'             => \$force);

if( $use_cache ) {
    &init_cache();
}
my $xs = XML::Simple->new;
my $csv = Text::CSV_XS->new ({ binary => 1, auto_diag => 1 });


while(<>) {
    chomp;
    s/^(\S+)\.(fq|fastq|fa|fasta)/$1/;
    my $query = $1;
    my $url = sprintf('esearch.fcgi?db=sra&tool=bioperl&term=%s',$query);
    my $output = get_web_cached($base,$url);
    warn("$output") if $debug;
    my $simplesum;
    eval {
	$simplesum = $xs->XMLin($output);
    };
    if( $@ ) {
	    delete_cache($base,$url);
	    next;
    }
    my $ids = $simplesum->{IdList}->{Id};
    if( ref($ids) !~ /ARRAY/ ) {
	$ids = [$ids];
    }
    for my $id ( @$ids ) {
	next if ! $id;
	if( $debug ) {
	    $url = sprintf('efetch.fcgi?retmode=text&db=sra&tool=bioperl&id=%s',$id);
	    $output = get_web_cached($base,$url);
	    eval { 
		$simplesum = $xs->XMLin($output);
	    };
	    warn(Dumper($simplesum)) if $debug;
	}
	$url = sprintf('elink.fcgi?dbfrom=sra&db=biosample&tool=bioperl&id=%s',$id);
	$output = get_web_cached($base,$url);
	eval { 
	    $simplesum = $xs->XMLin($output);
	};
	if( $@ ) {
	    delete_cache($base,$url);
	    next;
	}
	if( $simplesum->{LinkSet}->{LinkSetDb}->{Link} ) {
	    my $sampleref = $simplesum->{LinkSet}->{LinkSetDb}->{Link};
	    my $sampleid;
	    if( ref($sampleref) =~ /HASH/i ) {
		$sampleid = $sampleref->{Id};
	    } else {
		die("unsure how to handle a ",ref($sampleref), "for sample ID");
	    }
	    $url = sprintf('efetch.fcgi?db=biosample&tool=bioperl&id=%s',$sampleid);
	    $output = get_web_cached($base,$url);
	    my $biosample;
	    eval { 
		$biosample = $xs->XMLin($output);
	    };
	    if( $@ ) {
		delete_cache($base,$url);
		next;
	    }	    
	    warn("BioSAMPLE: ", Dumper($biosample)) if $debug;

	    my $ids = $biosample->{BioSample}->{Ids}->{Id};
	    if( ref($ids) !~ /ARRAY/ ) {
		$ids = [$ids];
	    }
	    my $strain_id;
	    for my $id ( @$ids ) {
		if( exists $id->{db_label} &&
		    $id->{db_label} eq 'Sample name' ) {
		    my $strain_id_str = $id->{content};
		    my @parsed = split(/_/,$strain_id_str);
		    # drop the last two _ separated bits
		    pop @parsed;
		    pop @parsed;
		    $strain_id = join("_", @parsed);
		}
	    }
	    if( $strain_id ) {
		print join("\t", $query, $strain_id),"\n";
	    } else {
		warn("No strain ID for $query\n");
	    }
	} else {
	    warn("no link for ", Dumper($simplesum),"\n");
	}
    }
    last if $debug;
}

sub init_cache {
    if( ! $cache_filehandle ) {
	mkdir($cache_dir) unless -d $cache_dir;
	$cache_filehandle = Cache::File->new( cache_root => $cache_dir);
    }
}

sub get_web_cached {
    my ($base,$url) = @_;
    if( ! defined $base || ! defined $url ) {
	die("need both the URL base and the URL stem to proceed\n");
    }
    unless( $use_cache ) {
	sleep $SLEEP_TIME;
	return get($base.$url);
    }
    my $val = $cache_filehandle->get($url);
    unless( $val ) {
	warn("$base$url not in cache\n") if $debug;
	$val = encode("utf8",get($base.$url));
	sleep $SLEEP_TIME;
	$cache_filehandle->set($url,$val,$cache_keep_time);
    }
    return decode("utf8",$val);
}

sub delete_cache {
    my ($base,$url) = @_;
    return unless $use_cache;

    if( ! defined $base || ! defined $url ) {
	die("need both the URL base and the URL stem to proceed\n");
    }

    $cache_filehandle->remove($base.$url);
}
