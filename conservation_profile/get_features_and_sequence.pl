#!/usr/bin/perl -w
use strict;

# Algorithm
# -- get coordinates of some features

use Bio::DB::SeqFeature::Store;
use Bio::SeqIO;
use Env qw(USER HOME);
use Getopt::Long;

my $Max_Intron = 800;

my ($user,$pass,$dbname,$host);
$host ='localhost';
my $odir = 'feature_alignments';
my $type = 'gene:NCBI_PASA';
my $debug = 0;
GetOptions(
	   'm|maxintron:i'=> \$Max_Intron,
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

my $transcriptout = Bio::SeqIO->new(-format => 'fasta',
				    -file   => ">$odir/$dbname.transcript.fasta");
my $utr3out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/$dbname.UTR3.fasta");
my $utr5out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/$dbname.UTR5.fasta");
my $cdsout = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/$dbname.CDS.fasta");
my $exonout = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/$dbname.exon.fasta");
my $intronout = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/$dbname.intron.fasta");

GENE: while( my $gene = $iter->next_seq ) {
    my ($gname) = $gene->get_tag_values('Note');
    if( defined $gname ) {
	$gname =~ s/^\"//;
        $gname =~ s/\"$//;
    }
    if( ! defined $gname || $gname !~ /^NCU/ ) {
	$gname = $gene->name;
    }
    warn("$gname ", $gene->location->to_FTstring(),"\n") if $debug;
    for my $mRNA ( $gene->get_SeqFeatures ) {
	warn(" $mRNA \n") if $debug;
	my ($transcript,$utr5,$utr3,$cds,$exon);
	my ($min,$max,$strand,$lastexon,$icount);	

	for my $exon ( sort { $a->start * $a->strand <=> 
			      $b->start * $b->strand } 
		       $mRNA->get_SeqFeatures('cds') ) {	    
	    warn("   exon $exon\n") if $debug;
	    $transcript .= $exon->seq->seq;
	    $min = $exon->start unless( defined $min && 
				       $min < $exon->start);
	    $max = $exon->end unless( defined $max && 
				     $max < $exon->end);
	    $strand = $exon->strand;

	    if( $lastexon ) {
		my ($intron_start, $intron_end);
		if( $strand > 0 ) { 
		    # last  current
		    # ------------>
		    # 1..5, 10..20
		    #  (1)    (2)
		    ($intron_start,$intron_end)= ( $lastexon->end+1,
						   $exon->start-1);
		} else {
		    # <------------
		    #  1..5   10..20
		    #  (2)     (1)
		    ($intron_start,$intron_end)= ( $exon->end +1,
						   $lastexon->start-1,
						   );
		}
		my $intronstr = sprintf($strand < 0 ? 
					"complement(%d..%d)" : "%d..%d",
					$intron_start,$intron_end);
		my $iseqobj = $dbh->segment($gene->seq_id,
					    $intron_start => $intron_end);
		if( ! defined $iseqobj ) {
		    warn("cannot find ",$gene->seq_id, ":$intron_start..$intron_end\n");
		    next GENE;
		} else {
		    $iseqobj = $iseqobj->seq;
		}
		
		my $iseq = $strand < 0 ? $iseqobj->revcom->seq : $iseqobj->seq;
		$intronout->write_seq(Bio::Seq->new
				      (-display_id=>sprintf("%s.i%d",
							    $mRNA->name,
							    ++$icount),
				       -desc => sprintf("%s:%s gene=%s",
							$gene->seq_id,
							$intronstr,
							$gname),
				       -seq => $iseq));
#		if( $iseqobj->length > $Max_Intron ) {
#		    warn("too long intron: $intronstr $gname ",$mRNA->name, ".$icount\n");
#		}
	    }
	    $lastexon = $exon;

	}
	if( ! defined $min || ! defined $max ) {
		warn("no min/max for ",$mRNA->name,"\n");
	}
	my $locstr = sprintf($strand < 0 ? "complement(%d..%d)" : "%d..%d",
			     $min,$max);
	
	$transcriptout->write_seq(Bio::Seq->new
				  (-display_id => 
				   sprintf("%s.cdna",
					   $mRNA->name),
				   -desc => sprintf("%s:%s gene=%s",
						    $gene->seq_id,
						    $locstr,
						    $gname),
				   -seq => $transcript));


	for my $utr (  sort { $a->start * $a->strand <=> 
			      $b->start * $b->strand } 
		       $mRNA->get_SeqFeatures('five_prime_utr') ) {
	    warn("   utr $utr\n") if $debug;
	    $utr5 .= $utr->seq->seq;
	    $min = $utr->start unless( defined $min && 
				       $min < $utr->start);
	    $max = $utr->end unless( defined $max && 
				     $max < $utr->end);
	    $strand = $utr->strand;
	}
	if( $utr5 ) {
	    $locstr = sprintf($strand < 0 ? "complement(%d..%d)" : "%d..%d",
			      $min,$max);
	    $utr5out->write_seq(Bio::Seq->new
				(-display_id => 
				 sprintf("%s.5utr",
					 $mRNA->name),
			     -desc => sprintf("%s:%s gene=%s",
					      $gene->seq_id,
					      $locstr,
					      $gname),
				 -seq => $utr5));
	}
	$min = undef;
	$max = undef;
	$strand = undef;
	for my $utr (  sort { $a->start * $a->strand <=> 
			      $b->start * $b->strand } 
		       $mRNA->get_SeqFeatures('three_prime_utr') ) {
	    warn("   utr $utr\n") if $debug;
	    $utr3 .= $utr->seq->seq;
	    $min = $utr->start unless( defined $min && 
				       $min < $utr->start);
	    $max = $utr->end unless( defined $max && 
				     $max < $utr->end);
	    $strand = $utr->strand;
	}
	if( $utr3 ) {
	    $locstr = sprintf($strand < 0 ? "complement(%d..%d)" : "%d..%d",
			      $min,$max);
	    $utr3out->write_seq(Bio::Seq->new
				(-display_id => 
				 sprintf("%s.3utr",
					 $mRNA->name),
				 -desc => sprintf("%s:%s gene=%s",
						  $gene->seq_id,
						  $locstr,
						  $gname),
				 -seq => $utr3));
	}
	$min = undef;
	$max = undef;
	$strand = undef;
	my $exonct = 1;
	for my $cdso (  sort { $a->start * $a->strand <=> 
			      $b->start * $b->strand } 
		       $mRNA->get_SeqFeatures('cds') ) {
	    warn("   cds $cdso\n") if $debug;
	    $cds .= $cdso->seq->seq;
	    $min = $cdso->start unless( defined $min && 
					$min < $cdso->start);
	    $max = $cdso->end unless( defined $max && 
				     $max < $cdso->end);
	    $strand = $cdso->strand;
	    $exonout->write_seq(Bio::Seq->new(-display_id => sprintf("%s.exon.%d",
								     $mRNA->name,
								     $exonct++,
								     ),
					      -desc => sprintf("%s:%s gene=%s",
							       $gene->seq_id,
							       $cdso->location->to_FTstring,
							       $gname),
					      -seq=>$cdso->seq->seq));
								     
	}
	if( $cds ) {
	    $locstr = sprintf($strand < 0 ? "complement(%d..%d)" : "%d..%d",
			      $min,$max);
	    $cdsout->write_seq(Bio::Seq->new
			       (-display_id => 
				sprintf("%s.cds",
					$mRNA->name),
				-desc => sprintf("%s:%s gene=%s",
						 $gene->seq_id,
						 $locstr,
						 $gname),
				-seq => $cds));
	}
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


