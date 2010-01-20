#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Env qw(USER HOME);
use Getopt::Long;

use lib "$HOME/src/bioperl/core";
use Bio::AlignIO;
# cutout all CDS in reference
#  - stitch together consensus sequence from all the individuals
#  - align in protein space
#  - project back to CDS, use MK-test with polarization


my $Datafile = '/projects/coccidiodies/variation/SnpReport_Cspp-3.txt';


my $Max_Intron = 800;

my ($user,$pass,$dbname,$host);
$dbname = 'gb_coccidioides_immitis_rs_2';
$host ='localhost';
my $odir = 'aln_CDS';
my $type = 'gene:CI2_FINAL_CALLGENES_3';
my $debug = 0;
my $refGenome = 'Ci_RS';
my $alndir = 'alignments/pecan_alignments';

my $set_missing_to_ref = 0;
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


open( my $fh => $Datafile ) || die "cannot open $Datafile: $!";
my $i = 0;
chomp(my $hdr = <$fh>);
$hdr =~ s/^\#//;
my %header = map { $_ => $i++ } split(/\t/,$hdr);

# A speedup to pre-index the by-species allele columns
my (@allele_cols,@sp_names, %sp_count);

for my $col ( grep { /^a_(\S+)/ } keys %header ) {
    unless ( $col =~ /^a_(C[IP])/ ) {
	warn("unknown allele column '$col'\n");
    }
    $sp_count{$1}++;
    $sp_count{'ALL'}++; # for total number of sample count, when data is not partitioned by species
    push @allele_cols, [ $1, $header{$col}, $col ];
}
# re-sort
@allele_cols = sort { $a->[0] cmp $b->[0] || 
		      $a->[2] cmp $b->[2] } @allele_cols;

my %snpByGene;
while(<$fh>) {
    chomp($_); # remove the trailing '\n'
    my @row = split(/\t/,$_); # split into an array
    my ($snpid,$lg,$start, $refallele,
	$alleles,$genes) = 
	    map { $row[ $header{$_} ] } qw(snpid linkageGroup pos 
					   refallele alleles genes);
    next unless $genes;
    $lg = 'Cimm_RS_2.'.$lg;
    for my $g ( split(/,/,$genes) ) {
	push @{$snpByGene{$genes}}, [ $snpid, $lg, $start, $refallele,
				      map { $row[$_->[1]] || '-' } @allele_cols];
    }
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );
my $iter = $dbh->get_seq_stream(-type => $type);


while(my $gene = $iter->next_seq ) {
    my $gname = $gene->name;
    warn("$gname\n");
    if( ! exists $snpByGene{$gname} ) {
	warn("skipping $gname, no SNPs\n");
    }
    my $written = 0;
    my %cds_seqs;
    my %locations;
    my $length;
    for my $mRNA ( $gene->get_SeqFeatures('mRNA') ) {
# what to do about ALT SPLICING?
	for my $exon ( sort { $a->start * $a->strand <=> 
				  $b->start * $b->strand } 
		       $mRNA->get_SeqFeatures('cds') ) {	    
	    warn("   exon $exon\n") if $debug;
	    my $exe = sprintf('sliceAlignment %s %s %s %d %d %s',$alndir, $refGenome, $exon->seq_id,
			      $exon->start-1,$exon->end, $exon->strand < 0 ? '-' : '+');
	    warn($exe,"\n") if $debug;
	    open(my $slicefh => "$exe |") || die $!;
	    my $alnin = Bio::AlignIO->new(-format => 'fasta',
					  -fh     => $slicefh);
	    my $ref_seq;
	    if( my $aln = $alnin->next_aln ) {
		for my $seq ( $aln->each_seq ) {
		    my $str = $seq->seq;
		    $str =~ s/\-//g;
		    $cds_seqs{$seq->id} .= $str;
		    $ref_seq = $str if( $seq->id eq $refGenome );
		    push @{$locations{$seq->id}}, $seq->description;
		    $length += $seq->length if $seq->id eq $refGenome;
		}
	    } else {
		next;
	    }
	    my @pop_seqs;
	    for my $col ( @allele_cols ) {
		push @pop_seqs, $ref_seq;
	    }
	    # now we will update the sequence
	    for my $snp ( @{$snpByGene{$gname}} ) {
	        my ($snpid, $lg, $start, $refallele,@alleles) = @$snp;
		# TODO need to check that the allele overlaps the exon!!
		for( my $i =0; $i < @alleles; $i++ ) {
		    my $allele = $alleles[$i];
		    unless( $set_missing_to_ref && $allele eq '-') { 
			# replace the allele
			my $offset = ( $exon->strand < 0 ) ? $exon->end : $exon->start;
			warn("replacing allele for ",$start, " for exon ",$exon->location->to_FTstring, " offset=$offset (",$offset-$start,") for $exon\n");
			warn("pop_seqs are $pop_seqs[$i]\n");
			substr($pop_seqs[$i],$offset-$start,1,$allele); 
		    }
		}
	    } 
	    my $i = 0;
	    for my $seq ( @pop_seqs ) {
		$cds_seqs{$allele_cols[$i]->[2]} .= $seq;
	    }
	}
	if( $length ) {
	    my $outc = Bio::SeqIO->new(-format => 'fasta',
				       -file   => ">$odir/$gname.CDS.fa");
	    my $outp = Bio::SeqIO->new(-format => 'fasta',
				       -file   => ">$odir/$gname.PEP.fa");
	    
	    for my $cds ( keys %cds_seqs ) {
		my $seq = Bio::PrimarySeq->new(-id  => $cds,
					       -seq => $cds_seqs{$cds},
					       -description => join(",", @{$locations{$cds} || []}));
		$outc->write_seq($seq);
		$outp->write_seq($seq->translate);
		$written = 1;
	    }
	}
	last; # probably want to only process 1st mRNA (random I guess)
    }
    last if $debug && $written;
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
