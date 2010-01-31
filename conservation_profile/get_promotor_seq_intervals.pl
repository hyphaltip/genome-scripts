#!/usr/bin/perl -w
use strict;

# Algorithm
# -- get coordinates of UTR features and 5' upstream in 1,2,3,4kb increments

use Bio::DB::SeqFeature::Store;
use Bio::SeqIO;
use Env qw(USER HOME);
use Getopt::Long;
use Bio::Location::Simple;
use Bio::SeqFeature::Generic;
use File::Temp qw(tempdir);
use constant SCALE => 1_000;


my ($user,$pass,$dbname,$host);
$host ='localhost';
my $odir = 'promotor_alignments';
my $type = 'gene:EMBL';
my $debug = 0;
my @increments = (1,2,3,4); # 1,2,3,4kb increments upstream
GetOptions(
	   'v|verbose!' => \$debug,
	   'u|user:s' => \$user,
	   'p|pass:s' => \$pass,
	   'host:s'   => \$host,
	   'db|dbname:s' => \$dbname,
	   'o|out|output:s' => \$odir,
	   't|type:s' => \$type,
	   );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

mkdir($odir) unless -d $odir;

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );

my $iter = $dbh->get_seq_stream(-type => $type);

my $tempdir = tempdir(CLEANUP => 1);
my $seen_intervals = Bio::DB::SeqFeature::Store->new(-adaptor => 'berkeleydb',
						     -dir     => $tempdir,
						     -create  => 1);
my %fh;
for my $size (@increments) {
    $fh{$size} =  Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$odir/$dbname.$type.$size.fasta");
}

GENE: while( my $gene = $iter->next_seq ) {
    my ($gname) = $gene->get_tag_values('locus_tag');
    if( defined $gname ) {
	$gname =~ s/^\"//;
        $gname =~ s/\"$//;
    }
    my $chrom_segment = $dbh->segment($gene->seq_id);
    my $seen_segment = $seen_intervals->segment($gene->seq_id);
    unless( defined $seen_segment ) {	
	$seen_intervals->store
	    (Bio::SeqFeature::Generic->new(-seq_id=>$gene->seq_id,
					   -strand => '+',
					   -start => 1,
					   -primary_tag  => 'scaffold',
					   -source_tag => 'chromosome',
					   -end   => $chrom_segment->length,
					   -tag   => {'Id' => $gene->seq_id,
						      'Name'=> $gene->seq_id}));
	$seen_segment = $seen_intervals->segment($gene->seq_id);
    }
    my $length = $chrom_segment->length;
    warn("$gname ", $gene->location->to_FTstring(),"\n") if $debug;
    my $start = $gene->start;
    my $end   = $gene->end;
    my $strand = $gene->strand;
    for my $size ( @increments ) {
	my $fh = $fh{$size};
	my @upstream;
	if( $strand < 0 ) {
	    last if $end == $length;
	    @upstream = ($end + SCALE,
			 $end + 1);	    
	    $upstream[0] = $length if( $upstream[0] > $length );
	} else {
	    last if $start == 1;
	    @upstream = ($start - SCALE,
			 $start - 1);
	    $upstream[0] = 1 if( $upstream[0] < 1 );
	}              
	
	my $segment = $dbh->segment($gene->seq_id,
				    sort { $a <=> $b } @upstream);
	
	if( $segment->length < 10 ) {
	    warn("skipping segment $segment, too short\n");
	}
	warn("segment is '$segment' for '", $gene->seq_id,"'\n") if $debug;
	my @features = $segment->features(-type => $type);
	
	if( @features ) {
	    warn("overlapping @features for @upstream of $segment for $gname\n");
	    last;
	}
	my $seqstr = $strand < 0 ? $segment->seq->revcom->seq :$segment->seq->seq;
	my $loc = $segment->location;
	$loc->seq_id($gene->seq_id);
	$loc->strand($strand);
	
	my @seenf = $seen_intervals->features(-seq_id =>$loc->seq_id, 
					      -start  => $loc->start, 
					      -end    => $loc->end,
					      -type   => 'PROMOTOR:SEEN');
	if(  @seenf ) {
	    warn("already covered interval: "
		,$seenf[0]->source_tag, " ",$seenf[0]->primary_tag, " ",
		$seenf[0]->location->to_FTstring(),"\n");
	    last;
	}		    
						 
	$fh->write_seq(Bio::PrimarySeq->new(-seq => $seqstr,
					    -display_id=>sprintf("%s.5UP_%dkb",
								 $gname,
								 $size),
					    -desc => sprintf("%s:%s gene=%s",
							     $gene->seq_id,
							     $loc->to_FTstring(),
							     $gname)
					    ));
	my $f = Bio::SeqFeature::Generic->new(-seq_id   => $gene->seq_id,
					   -source_tag   => 'SEEN',
					   -primary_tag     => 'PROMOTOR',
					   -location => $loc);
	

	$seen_intervals->store($f);

	if( $strand > 0 ){ 
	    $start = $upstream[0];	    
	} else { 
	    $end = $upstream[0];	    
	}
	last if( ($upstream[0] == 1 && $strand > 1) || 
		 ($upstream[1] > $length && $strand < 0));
    }
	
    last if $debug;
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


