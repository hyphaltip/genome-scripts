#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Env qw(HOME);
my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $src = 'gene';
my $output;

GetOptions(
	   'v|verbose!'  => \$debug,
	   'u|user:s'    => \$user,
	   'p|pass:s'    => \$pass,
	   'host:s'      => \$host,
	   'db|dbname:s' => \$dbname,
	   
	   's|src:s'     => \$src,
	   'o|output:s'  => \$output,
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

$output ||= sprintf("%s-%s.intron.fa",$dbname,$src);
my $ofh = Bio::SeqIO->new(-format => 'fasta',
			  -file   => ">$output");
my @names;
my $iter = $dbh->get_seq_stream(-type => 'scaffold:chromosome');
while( my $segment = $iter->next_seq ) {
#    my $segment = $dbh->segment($seg->seq_id);
    my $count = 0;
    warn("seg is $segment\n");
    # FIX ME HERE
    for my $gene ($segment->get_SeqFeatures(-type => 'gene') ) { #(-type => "$src:") ) {
	warn("gene is $gene\n");
	my $gname = $gene->name;
	my @mrna = $gene->get_SeqFeatures('mRNA');
	for my $mRNA ( @mrna ) {
	    my ($id) = $mRNA->load_id;
	    my @exons = $mRNA->get_SeqFeatures('cds');
	    my $lastexon;
	    my $i = 1;
	    for my $exon ( sort { ($a->start*$a->strand) <=> 
				      ($b->end*$b->strand) } 
			   @exons ) {
		if( $lastexon ) {
		    my ($seqid,$s,$e) = ($exon->seq_id);
		    my $loc;
		    if( $exon->strand > 0 ) {			
			($s,$e) = ($lastexon->end+1, $exon->start-1);
			$loc = Bio::Location::Simple->new(-start => $s,
							  -end   => $e);
		    } else {
			($s,$e) = ($lastexon->start-1,$exon->end+1);	
			$loc = Bio::Location::Simple->new(-start => $e,
							  -end   => $s,
							  -strand => $exon->strand
			    );    
		    }
		    my $iseq = $dbh->segment($seqid,$s,$e);
		    if( ! $iseq ) {
			warn("cannot get seq for $seqid\n");
		    } else {
			$iseq = $iseq->seq->seq;
		    }

		    my $intron = Bio::PrimarySeq->new
			(-seq => $iseq,
			 -id => "$prefix:$id.i".$i++,
			 -desc => sprintf("%s %s:%s",
					  $gname,
					  $seqid,$loc->to_FTstring));
		    if( $intron->length < 4) {
			warn("intron ",$intron->length, " too short for: ", 
			     $intron->display_id,"\n");
		    } else {
			$ofh->write_seq($intron);
		    }
		}
		$lastexon = $exon;
		last if $debug && $count++ > 10;
	    }
	}
    }
    last;
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
