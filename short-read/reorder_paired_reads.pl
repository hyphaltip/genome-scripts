#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;
use File::Temp qw(tempdir);
my $tempbase;
my ($min_length, $min_qual,$min_percent) = (28,30,70);
my $offset = 33; # Sanger offset is 33, Illumina is 64
my $outdir = ".";
my $debug = 0;
GetOptions(
	   'v|verbose!'  => \$debug,
	   'l|length:i'  => \$min_length,
	   'q|qual:i'    => \$min_qual,
	   'p|percent:i' => \$min_percent,
	   't|tempbase:s'=> \$tempbase,
	   'o|outdir:s'  => \$outdir,
	   'f|offset:i'  => \$offset,
	   );

# don't use this on files with more than 1M reads it will be really slow
# split the files with fastq_split

my ($ext) = "(fq|fastq|txt|seq)";

my %sets;
for my $file (@ARGV) {    
    $file = File::Spec->rel2abs($file);
    my ($vol,$dir,$fname) = File::Spec->splitpath($file);
    my @name = split(/\./,$fname);
    my ($ext) = pop @name;
    my $f = join(".",@name);
    if( $f =~ /(\S+_s_\d+)\_(\d+)\.p(\d+)$/ ||
	$f =~ /(\S+)\_(\d+)\.p(\d+)$/ ) {
	$sets{"$1.$3"}->{$2} = $file;
    } elsif( $f =~ /(\S+)\.p(\d+)$/ ) {
	$sets{$1}->{$2} = $file;
    }
}
my %conditions;
for my $set ( keys %sets ) {
    my ($base) = split(/\./,$set);
    my $temp = tempdir(DIR => $tempbase, CLEANUP => 1);
    warn "temp is $temp\n" if $debug;
    for my $t ( keys %{$sets{$set}} ) {	
	warn "$set $t ",$sets{$set}->{$t},"\n" if $debug;
	my $exe = sprintf(<<EOL
fastq_quality_filter -v -Q %d -q %d -p %d -i %s | fastq_quality_trimmer -v -Q %d -t %d -l %d -o %s
EOL
,
			  $offset, $min_qual, $min_percent, $sets{$set}->{$t},
			  $offset, $min_qual, $min_length, 
			  "$temp/$set.$t\_\_$base");	
	`$exe`;
    }
    opendir(DIR,$temp);
    my %tset;
    for my $outf ( readdir(DIR) ) {
	next if $outf =~ /^\./;
	if( $outf =~ /(\S+)\.(\d+)__(\S+)/ ){	    
	    my ($base,$lane,$cond) = ( $1,$2,$3 );
	    warn "$base -> $outf\n" if $debug;
	    $tset{$cond}->{$lane} = [$base,"$temp/$outf"];
	}	
    }
    for my $cond ( keys %tset ) {
	my $ofh_cond_1 = $conditions{"$cond.1"};	
	if( ! defined $ofh_cond_1 ) {
	    open($ofh_cond_1 => ">$outdir/$cond\_1.fq") || die $!;
	    $conditions{"$cond.1"} = $ofh_cond_1;
	}

	my $ofh_cond_2 = $conditions{"$cond.2"};	
	if( ! defined $ofh_cond_2 ) {
	    open($ofh_cond_2 => ">$outdir/$cond\_2.fq") || die $!;
	    $conditions{"$cond.2"} = $ofh_cond_2;
	}

	my $ofh_cond_unpair = $conditions{"$cond.unpair"};	
	if( ! defined $ofh_cond_unpair ) {
	    open($ofh_cond_unpair => ">$outdir/$cond\_unpair.fq") || die $!;
	    $conditions{"$cond.unpair"} = $ofh_cond_unpair;
	}
	my %alldata;

	for my $lane ( keys %{$tset{$cond}} ) {
	    my ($b,$fname_filtered) = @{$tset{$cond}->{$lane}};
	    warn "$cond $lane base=$b file=$fname_filtered\n";
	    open(my $infh => $fname_filtered) || die $!;
	    while(<$infh>) {
		my $record = $_;
		my $name;
		if( /^@(\S+)\#\d+\/([12])/ ||
		    /^@(\S+\.\d+)\s+\S+:\d+\/([12])/  ) {
		    $name = $1;
		    if( $2 != $lane ) {
			chomp($_);
			warn("lane ($lane) does not match read name ($_)\n");
		    }
		} elsif( /^@(\S+):[YN]/) {
   	 		$name = $1;	
		} else{
			warn("cannot parse name out for $_");
			next;
		}
		for ( 0..2) {
		    $record .= <$infh>;
		}
		if( ! defined $name ) {
		    warn("no name for $record\n");
		    die;
		} elsif( ! defined $lane ) {
		    warn("no lane for $record\n");
		    die;
		}
		$alldata{$name}->{$lane} = $record;
	    }
	}
	for my $readname ( sort keys %alldata ) {
	    my $d = $alldata{$readname};
	    if( exists $d->{1} && exists $d->{2} ) {
		print $ofh_cond_1 $d->{1};
		print $ofh_cond_2 $d->{2};
	    } elsif( exists $d->{1} ) {
		print $ofh_cond_unpair $d->{1};
	    } elsif( exists $d->{2} ) {
		print $ofh_cond_unpair $d->{2};
	    }
	}
    }
}


