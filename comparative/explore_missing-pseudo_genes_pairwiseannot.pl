#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::DB::Fasta;


my $ref_protein_db; # protein DB of the proteins from the ref genome
my $ref_cds_db; # protein DB of the proteins from the ref genome
my $target_genome;  # genomic (FASTA) of the target DB

my $outdir = 'pseudo_search';
my $infile;
my $debug;

GetOptions('r|ref:s'      => \$ref_protein_db,
	   'c|cds:s'      => \$ref_cds_db,
	   't|target:s'   => \$target_genome,
	   'o|out|outdir:s'=> \$outdir,
	   'i|in|infile:s' => \$infile,
	   'v|verbose:s'   => \$debug,
	   );


my $refdb = Bio::DB::Fasta->new($ref_protein_db);
my $cdsdb = Bio::DB::Fasta->new($ref_cds_db);
my $tdb = Bio::DB::Fasta->new($target_genome);


$infile ||= shift @ARGV;
my ($stem);
if( $infile =~ /-vs-(\S+)\.tab/ ){
    $stem = $1;
} else {
    $stem = $infile;
    $stem =~ s/\.(\S+)$//;
}

mkdir($outdir);
mkdir("$outdir/$stem");

open(my $fh => $infile) || die "cannot open $infile: $!";
my $header = <$fh>;
my $i = 0;
my %header = map { $_ => $i++ } split(/\s+/,$header);
my @data;
my @targets;


# DATA LOOKS LIKE THIS
# GENE_FROM	MRNA_FROM	CHROM_FROM	START_FROM	STOP_FROM	STRAND_FROM	GENE_TO	MRNA_TO	START_TO	STOP_TO	STRAND_TO	SINGLE_MATCH
# NCU09290 NCU09290T0 supercont10.7 4192751 4194673 -1 maker-assembled104-augustus-gene-0.62 maker-assembled104-augustus-gene-0.62-mRNA-1 assembled104 43470 45152 -1 yes
# NCU09289 NCU09289T0 supercont10.7 4191519 4192340 1 NO_GENES_IN_INTERVAL
#  

while(<$fh>) {
    my @row = split;    
    if( $row[ $header{GENE_TO} ] eq 'NO_GENES_IN_INTERVAL' ) {
	push @targets, scalar @data; # index of current row
    }
    push @data, [@row];
}

for my $target ( @targets ) {
    next if $target == 0; # skip beginning, want to be flanked by two genes
    my ($before,$here,$after) = map { $data[$_] } ($target-1,$target, $target+1);
    if( ! defined $before->[ $header{MRNA_TO} ] ||
	! defined $after->[ $header{MRNA_TO} ] ) {
	# skip if there isn't a gene before or after
	# this isn't really the coumn but it will have this data in that next col
	next;
    }
    if( $before->[ $header{CHROM_TO} ] ne $after->[ $header{CHROM_TO}] ||
	$before->[ $header{CHROM_FROM} ] ne $here->[$header{CHROM_FROM} ] || 
	$after->[ $header{CHROM_FROM} ] ne $here->[$header{CHROM_FROM} ] ) {
	# skip if contigs change
	next;
    }
    my $gene_from = $here->[ $header{GENE_FROM} ];
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$outdir/$stem/$gene_from.genomic.fa");
    
#    print "<< ",join(" ", @$before),"\n","  ",join(" ",@$here),"\n",">> ",join(" ",@$after),"\n";
    my ($region_start, $region_end )  = ($before->[ $header{STOP_TO} ]+1, 
					 $after->[ $header{START_TO} ] -1);
    if( $here->[ $header{STRAND_FROM} ] == -1 ) {
	# flip-flop so we get rev-strand genomic
	($region_end,$region_start) = ($region_start,$region_end);
    }
    my $seqstr = $tdb->seq($before->[ $header{CHROM_TO} ],$region_start,$region_end );
    my $seq = Bio::Seq->new(-id => sprintf("%s_%s",
					   $here->[ $header{GENE_FROM} ],
					   $stem),
			    -seq => $seqstr,
			    -description => sprintf("%s:%d..%d",
						    $before->[ $header{CHROM_TO} ],
						    $region_start,$region_end));
    $out->write_seq($seq);
    $out = Bio::SeqIO->new(-format => 'fasta',
			   -file   => ">$outdir/$stem/$gene_from.protein.fa");
    my $pepseq = $refdb->get_Seq_by_acc($here->[$header{MRNA_FROM} ]);
    if( ! defined $pepseq ) { 
	warn("cannot find $pepseq\n");
    }
    $out->write_seq($pepseq);
    $out = Bio::SeqIO->new(-format => 'fasta',
			   -file   => ">$outdir/$stem/$gene_from.cds.fa");
    my $cdsseq = $cdsdb->get_Seq_by_acc($here->[$header{MRNA_FROM} ]);
    if( ! defined $cdsseq ) { 
	warn("cannot find $cdsseq\n");
    }
    $out->write_seq($cdsseq);
    last if $debug;
}
