#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max shuffle);
use Bio::DB::SeqFeature;
use Bio::Seq;

use Env qw(USER HOME);
my %uncompress = ( 'gz' => 'zcat',
		   'bz2'=> 'bzcat',
		   );
my $debug;
my $gff_dir;
my ($dbname,$transcript_source);
my $dsn = 'dbi:mysql';
my $adaptor = "DBI::mysql";
my ($input,$output);
my $variant_source;
my $flatfile_adaptor = "berkeleydb";
my $min_coverage = 5;
GetOptions(
	   'debug|verbose!'     => \$debug,

	   'h|help'      => sub {  # help menu 
               exec('perldoc', $0);
               exit(0);
           },

	   'f|gff:s'        => \$gff_dir,
	   # or
	   'db|dbname:s'    => \$dbname,
	   'dsn:s'          => \$dsn,
	   'adaptor:s'      => \$adaptor,
	   's|source:s'        => \$transcript_source,
	   'vs|varsource:s'  => \$variant_source,
	   'i|input=s'      => \$input,

	   'm|min|coverage:i'=> \$min_coverage,
	   'o|output:s'      => \$output,
	   );

unless( $input ) {
    die( "must provide an input pileup file\n");
}
if( ! defined $output ) {
    $output = "$input.gff3";
}

unless( $variant_source ) {
    # if not specified
    # take first bit before first '.' as the variant_source
    ($variant_source) = split(/\./,$input);
}
my $fh;
if( $input =~ /\.(gz|bz2)$/ ) {
    open($fh => "$uncompress{$1} $input |") || die "cannot open $input with uncompress on $1\n";
} else {
    open($fh => "< $input") || die "cannot open $input: $!";
}

my $dbh;
my @group_fields = qw(ID Reference_seq Variant_seq 
		      Variant_reads Genotype Total_reads Variant_effect); 
my @src;
if( $dbname  && $transcript_source ) {
    $dbh = Bio::DB::SeqFeature::Store->new(-adaptor=> $adaptor,
					  -dsn    => "$dsn:$dbname",
					  );    
    @src= ("mRNA:$transcript_source");
} elsif( $gff_dir ) {
    $dbh = Bio::DB::SeqFeature::Store->new(-adaptor=> $flatfile_adaptor,
					  -dir    => $gff_dir,
					   );
    @src = ("mRNA:$transcript_source" || 'gene');
} else {
    die("Must provide either dbname and gene source field for Bio::DB::SeqFeature database OR gff file (and genome fasta)\n");
}
open(my $ofh => ">$output") || die "cannot open $output for writing\n";
warn("source is @src\n") if $debug;
while(<$fh>) {
    chomp;
    my ($chrom,$pos,$ref_base,$read_allele,
	$consensus_qual, $snp_qual,
	$RMS_map_qual, $num_reads,@row) = split(/\t/,$_);
    my %group = ('ID' => sprintf("%s_%d",$chrom,$pos),
		 'Reference_seq' => $ref_base,
		 'Variant_seq' => $read_allele,
		 'Total_reads' => $num_reads);
    
    my $segment = $dbh->segment($chrom, $pos-1 => $pos+1);
    my @features = $segment->features(-type => [@src]);
    
    my $var_length = 1; # how many bases
    my $var_type;
    my $type;
    if( $ref_base eq '*' ) {	# indel
	my ($first_allele, $second_allele, $num_reads_sup_first, 
	    $num_reads_sup_second, 
	    $num_reads_containg_indels) = @row;	

	next unless $num_reads_sup_second > $min_coverage;
#	print "@features overlap $chrom:$pos:$ref_base --> $snp_qual, $consensus_qual $read_allele $num_covering_reads($num_reads_sup_first:$num_reads_sup_second)\n";
	$group{'Variant_reads'} = $num_reads_sup_second;
	
	for my $r ( split(/\//,$read_allele )) {
	    next if $r eq '*';
	    my $indeldir = substr($r,0,1,''); # get the - or + and strip it out
	    if( $indeldir eq '-' ) {
		if( defined $var_type && $var_type ne 'deletion' ) {
		    $var_type = 'complex_substitution';
		    last;
		}
		$var_type = 'deletion';
	    } elsif( $indeldir eq '+' ) {
		if( defined $var_type && $var_type ne 'insertion' ) {
		    $var_type = 'complex_substitution';
		    last;
		}
		$var_type = 'insertion';
	    } else {
		warn("unknown indeldir for '$r' - $indeldir ($read_allele)\n");
	    }
	    $var_length = max($var_length,length($r)); # longest is the event len
	}
	$type = 'nucleotide_'.$var_type;
    } else {
	my ( $num_read_bases, $read_qual, $aln_qual) = @row;
	$type = 'SNV';
	$var_type = 'point_mutation';
    }
    if( @features ) {
	
	# determine if this is AA changing
	for my $mRNA ( @features ) {
	    my ($point) = $dbh->segment($chrom,$pos,$pos+1);	    
	    # print "mRNA is ",$mRNA->location->to_FTstring, " pos is $pos\n";
	    my ($cds,$altcds);
	    my $name = $mRNA->name;
	    
	    my ($CDS_Overlap) = $point->features(-type => ["CDS:$transcript_source"]);
	    my $cds_offset;
	    if( ! $CDS_Overlap ) {
		# check for splice-site
		$group{'Variant_effect'}= sprintf("%s %d %s %s",
						  $var_type, $var_length,
						  'intron',
						  $name);
		last;
	    } elsif( $type ne 'SNV' ) {
		# this is an indel - it must be a frameshift
		$group{'Variant_effect'} = sprintf("%s %d %s %s",
						   'frameshift_sequence_variation',
						   $var_length,
						   $mRNA->primary_tag,
						   $name);
		last;				
	    }
	    
	    my $n =1;
	    for my $CDS ( sort { $a->start * $a->strand <=> 
				     $b->start * $b->strand } 
			  $mRNA->get_SeqFeatures('CDS') ) {
		my $this_CDS = $dbh->seq($chrom,$CDS->start => $CDS->end);
		my $CDS_seq = Bio::Seq->new(-seq => $this_CDS);
		$cds .= $CDS->strand > 0 ? $this_CDS : $CDS_seq->revcom->seq;
		# change it to the variant allele		
		if($CDS->load_id eq $CDS_Overlap->load_id ) {
		    warn($CDS->start,"..",$CDS->end, " for ", $CDS->location->to_FTstring(), "\n") if $debug;
		    my $cds_offset = $pos - $CDS->start;
		    substr($this_CDS,$cds_offset,1,$read_allele);
		    warn("this_cds - ($n) $cds_offset for $pos in $this_CDS ($CDS == $CDS_Overlap)\n") if $debug;
		    
		}
		$CDS_seq = Bio::Seq->new(-seq => $this_CDS);
		$altcds .= $CDS->strand > 0 ? $this_CDS : $CDS_seq->revcom->seq;
		$n++;
	    }
	    warn("$ref_base --> $read_allele\n") if $debug;

	    my $pep = Bio::Seq->new(-seq=>$cds)->translate->seq;
	    my $altpep = Bio::Seq->new(-seq=>$altcds)->translate->seq;
	    
	    if( $debug ) {
		open(my $t1 => ">test.fa") || die $!;
		print $t1 ">orig_$name\n$cds\n";
		print $t1 ">alt\n$altcds\n";
		open(my $t2 => ">test2.pep") || die $!;
		print $t2 ">orig_$name\n$pep\n";
		print $t2 ">alt\n$altpep\n";
	    }
	    my $effect;
	    $pep    =~ s/\*$//;
	    $altpep =~ s/\*$//;
	    if( $altpep ne $pep ) {		
		if( $altpep =~ /\*/ ) {
		    $effect = 'sequence_variant_causing_nonsense_codon_change_in_transcript';
		} else {
		    $effect = 'sequence_variant_causing_missense_codon_change_in_transcript';
		}
	    } else {
		 $effect = 'sequence_variant_causing_synonymous_change_in_transcript';
	    }
	   
	    $group{'Variant_effect'} = sprintf("%s %d %s %s",
					       $effect,
					       $var_length,
					       $mRNA->primary_tag,
					       $name);
#	    exit;# if $CDS_Overlap->strand > 0;
	}
    } else {
	$group{'Variant_effect'}= sprintf("%s %d %s",
					  $var_type, $var_length,
					  'intergenic_region');
    }
# strand is '+' I guess
# frame is '.' I guess

    
    my $group = join(";", grep { defined }
		     map { exists $group{$_} ? sprintf("%s=%s",$_,$group{$_}) : undef } @group_fields);
    print $ofh join("\t", $chrom,$variant_source,$type,
	       $pos, $pos+1, $snp_qual, '+',
	       '.',$group), "\n";

}

__END__

=head1 NAME

pileup2gvf - samtools pileup to GVF (genome variant format)

=head1 SYNOPSIS

 pileup2gvf.pl --gff annotation.gff3 --genome genome.fasta -i project.pileup

 pileup2gvf.pl --gff annotation.gff3 --genome genome.fasta -i project.pileup.gz

 pileup2gvf.pl --db projectgenome -s GenBank -i project.pileup.gz

=head1 DESCRIPTION

 -i --input   - [r] the pileup input. Compressed (.gz or .bz2) files are automatically uncompressed on the fly, 
                requiring zcat or bzcat to be installed. 
 -o --output  - the output file to store this in.

B<Annotation data source>

One of these sets is required to describe where the annotation is located. Either the Bio::DB::SeqFeature database

  -db --dbname - the relational database name (assumed dsn prefix is dbi:mysql, can be changed with --dsn option) 
  -s --source  - the source field from the GFF3 that was loaded

Or the GFF3 file
  -f --gff     - the GFF3 file with the annotation (gene/mRNA/CDS/exon)
  -g --genome --dna - the genome fasta (if not encoded in GFF3 file
 
=head1 DESCRIPTION

SAMtools pileup format is documented here:
http://samtools.sourceforge.net/cns0.shtml

GVF is documented in part in the recent paper:
 http://dx.doi.org/10.1186/gb-2010-11-8-r88

And SO (Sequence Ontology) is documented here http://www.sequenceontology.org/

Also see the mutations event calling page for more info on event types.
http://www.ebi.ac.uk/mutations/recommendations/mutevent.html

=head1 AUTHOR

Jason Stajich - jason.stajich [at] ucr.edu or jason [at] bioperl.org

=cut
