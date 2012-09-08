#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use List::Util qw(sum);
use Env qw(HOME);
my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $max_5_offset = 6;
my $max_3_offset = 3;
my $src = 'gene:NC10_CALLGENES_FINAL_5';
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
my $src_sanitize = $src;
$src_sanitize =~ s/:/_/g;
$output ||= sprintf("%s-%s.intron.fa",$dbname,$src_sanitize);

my $table = sprintf("%s-%s.intron.tab",$dbname,$src_sanitize);
my $ofh = Bio::SeqIO->new(-format => 'fasta',
			  -file   => ">$output");
open(my $ofhTable => ">$table") || die "cannot create $table: $!";
print $ofhTable join("\t", qw(gene mRNA intron 5SS 5SS_location 
                              BP BP_location 
                              3SS 3SS_location 
                              bpA_to_3ss intron_length)),"\n";

my $freqTb = sprintf("%s-%s.intron_SSfreq.tab",$dbname,$src_sanitize);
open(my $freqTable => ">$freqTb") || die "cannot create $freqTb: $!";

my @names;
my $iter = $dbh->get_seq_stream(-type => $src);
my %freq;
my $count = 0;
my $icount = 0;
while( my $gene = $iter->next_seq ) {
    my $gname = $gene->name || $gene->load_id;
    warn("gene is $gene name is $gname\n") if $debug;
    my @mrna = $gene->get_SeqFeatures('mRNA');
    for my $mRNA ( @mrna ) {
	my ($mRNA_id) = $mRNA->name || $mRNA->load_id;
	warn("mRNA is $mRNA name is $mRNA_id\n") if $debug;
	my @exons = $mRNA->get_SeqFeatures('cds');
	if( ! @exons ) {
	    @exons = $mRNA->get_SeqFeatures('exon');
	}
	warn("exon count is ",scalar @exons, "\n") if $debug;
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
						      -end   => $e,
						      -strand => 1,
			);
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
		     -id => "$mRNA_id.i".$i,
		     -desc => sprintf("%s %s:%s",
				      $gname,
				      $seqid,$loc->to_FTstring));
		if( $intron->length < 4) {
		    warn("intron ",$intron->length, " too short for: ", 
			 $intron->display_id,"\n");
		} else {
		    $ofh->write_seq($intron);
		    for my $offset ( 0..($max_5_offset-1) ) {
			$freq{'5SS'}->[$offset]->{substr($iseq,$offset,1)}++;
		    }
		    # -3,-2,-1
		    for my $offset ( map { -1 * $_ } 
				     1..$max_3_offset ) {
			$freq{'3SS'}->[$max_3_offset + 
				       $offset]->{substr($iseq,$offset,1)}++;
		    }
		    print $ofhTable join("\t", $gname, $mRNA_id,$i,
					 substr($iseq,0,$max_5_offset),
					 $loc->strand > 0 ? $loc->start : $loc->end,
					 '','', # Branchpoint
					 substr($iseq,-1 * $max_3_offset,$max_3_offset),
					 $loc->strand > 0 ? $loc->end : $loc->start,
					 '', #bpAdistanc to 3'ss
					 $loc->length),"\n";
		    $icount++;
		}
		    $i++;
	    }
	    $lastexon = $exon;
	}
    }
    last if $debug && $count++ > 10;
}
my @bases = qw ( G A T C);

for my $f ( sort keys %freq ) {
    print $freqTable "$f:\n";    
    my $max = scalar @{$freq{$f}} - 1;
    for my $letter ( @bases ) {
	print $freqTable join("\t", $letter, map { $freq{$f}->[$_]->{$letter} || 0 } 0..$max),"\n";
    }
    
    print $freqTable "\n";
    for my $letter ( @bases ) {
	print $freqTable join("\t", $letter, 
			      map { sprintf("%.2f",
					    ($freq{$f}->[$_]->{$letter} || 0) / 
					    $icount ) } 0..$max),"\n";
    }
    print $freqTable "\n";
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
