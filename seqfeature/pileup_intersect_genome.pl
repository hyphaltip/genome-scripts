#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Env qw(HOME);
my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my @src = 'gene:NCBI_PASA';
my ($input,$output);
my $min_coverage = 3;
GetOptions(
	   'v|verbose!' => \$debug,
	   'u|user:s' => \$user,
	   'p|pass:s' => \$pass,
	   'host:s'   => \$host,
	   'db|dbname:s' => \$dbname,
	   
	   's|src:s'     => \@src,
	   'i|input:s'   => \$input,
	   'o|output:s'  => \$output,
	   );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

$input ||= shift @ARGV;

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );

if( $output && $output ne '-' ) { 
    open($output => ">$output" ) || die $!;
} else {
    $output = \*STDOUT;
}

open(my $fh => $input) || die $!; 

# assuming samtools pileup 

# Print the alignment in the pileup format. In the pileup format, each
# line represents a genomic position, consisting of chromosome name,
# coordinate, reference base, read bases, read qualities and
# alignment mapping qualities. Information on match, mismatch, indel,
# strand, mapping quality and start and end of a read are all encoded
# at the read base column. At this column, a dot stands for a match to
# the reference base on the forward strand, a comma for a match on the
# reverse strand, ‘ACGTN’ for a mismatch on the forward strand and
# ‘acgtn’ for a mismatch on the reverse strand.  A pattern
# ‘\+[0-9]+[ACGTNacgtn]+’ indicates there is an insertion between this
# reference position and the next reference position. The length of
# the insertion is given by the integer in the pattern, followed by
# the inserted sequence. Similarly, a pattern ‘-[0-9]+[ACGTNacgtn]+’
# represents a deletion from the reference.  The deleted bases will be
# presented as ‘*’ in the following lines. Also at the read base
# column, a symbol ‘^’ marks the start of a read segment which is a
# contiguous subse- quence on the read separated by ‘N/S/H’ CIGAR
# operations. The ASCII of the character following ‘^’ minus 33 gives
# the mapping quality. A symbol ‘$’ marks the end of a read segment.

                 
# If option -c is applied, the consensus base, Phred-scaled consensus
# quality, SNP quality (i.e. the Phred-scaled probability of the
# consensus being identical to the reference) and root mean square
# (RMS) mapping quality of the reads covering the site will be
# inserted between the ‘reference base’ and the ‘read bases’
# columns. An indel occupies an additional line. Each indel line
# consists of chromosome name, coordinate, a star, the genotype,
# consensus quality, SNP quality, RMS mapping quality, # covering
# reads, the first alllele, the second allele, # reads supporting the
# first allele, # reads supporting the second allele and # reads
# containing indels different from the top two alleles.

while(<$fh>) {
    chomp;
    my ($chrom,$pos,$ref_base,@row) = split(/\t/,$_);
    my $segment = $dbh->segment($chrom, $pos-1 => $pos+1);
    my @features = $segment->features(-type => [@src]);
    if( $ref_base eq '*' ) { # indel
	my ($read_allele, $consensus_qual, $snp_qual, $RMS_mapping_qual,
	 $num_covering_reads, $first_allele, $second_allele, $num_reads_sup_first, $num_reads_sup_second, 
	 $num_reads_containg_indels) = @row;

	next unless $num_reads_sup_second > $min_coverage;

	print "@features overlap $chrom:$pos:$ref_base --> $snp_qual, $consensus_qual $read_allele $num_covering_reads($num_reads_sup_first:$num_reads_sup_second)\n";
	
    } else {
	my ( $num_read_bases, $read_qual,$aln_qual) = @row;
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
