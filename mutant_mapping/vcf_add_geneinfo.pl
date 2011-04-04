#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw(min max shuffle);
use Bio::DB::SeqFeature;
use Bio::Seq;
use Bio::Tools::IUPAC;
use Bio::Coordinate::GeneMapper;


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
	   's|source:s'     => \$transcript_source,
	   'vs|varsource:s' => \$variant_source,
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
		      Variant_reads Genotype Total_reads Variant_effect Codon Ref_Codon Variant_Codon Gene_strand
		      ); 
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
my $row = 0;
my (%sample, %col, %cnt);
while(<$fh>) {
    next if /^\##/;
    if (/^#/) {
      my @t = split("\t");
      for my $i (9 .. $#t) {
	my $k = $sample{$t[$i]};
	if ($k) {
	  push(@{$col{$k}}, $i);
	  push(@{$cnt{$k}}, 0, 0);
	}
      }
    } else {
      chomp;
      my ($chrom,$pos,$id,$ref_base,$alt, $qual,$filter, $info, $format, @genotypes) = split(/\t/,$_);

      next unless ( $alt ne '.' && $alt !~ /,/ );
      $row++;
      my @read_alleles;
      my %group = ('ID' => sprintf("%s_%d",$chrom,$pos),
		   'Reference_seq' => $ref_base,
		   'Variant_seq' => $alt,
		   'Quality' => $qual);
      my $segment = $dbh->segment($chrom, $pos-1 => $pos+1);
      my @features = $segment->features(-type => [@src]);
      my $var_length = 1;	# how many bases
      my $var_type;
      my $type;
      #my ( $num_read_bases, $read_qual, $aln_qual) = @row;
      $type = 'SNV';
      my $ambiseq = Bio::Seq->new(-seq => $alt, -alphabet => 'dna');	
      my $stream  = Bio::Tools::IUPAC->new(-seq => $ambiseq);
      while (my $uniqueseq = $stream->next_seq()) {
	my $char = $uniqueseq->seq;	    
	push @read_alleles, $char unless $char eq $ref_base;
      }
      $group{'Variant_seq'} = join(',', @read_alleles);
      $var_type = 'point_mutation';
      if ( @features ) {
	# determine if this is AA changing
	for my $mRNA ( sort {$b->length <=> $a->length} @features ) {
	  $group{'Gene_strand'} = $mRNA->strand;
	  my ($point) = $dbh->segment($chrom,$pos,$pos);	    
	  # print "mRNA is ",$mRNA->location->to_FTstring, " pos is $pos\n";

	  my $name = $mRNA->name;
	    
	  my ($CDS_Overlap) = $point->features(-type => ["CDS:$transcript_source"]);
	  my $cds_offset = 0;
	  my $cds = '';
	  if ( ! $CDS_Overlap ) {
	    # check for splice-site
	    $group{'Variant_effect'}= sprintf("%s %d %s %s",
					      $var_type, $var_length,
					      'intron',
					      $name);
	    last;
	  } elsif ( $type ne 'SNV' ) {
	    # this is an indel - it must be a frameshift
	    $group{'Variant_effect'} = sprintf("%s %d %s %s",
					       'frameshift_sequence_variation',
					       $var_length,
					       $mRNA->primary_tag,
					       $name);
	    last;
	  }
	    
	  my $gene_coords = Bio::Coordinate::GeneMapper->new(-in => 'chr',
							     -out=> 'cds',
							     -exons => [$mRNA->get_SeqFeatures('CDS')],
							    );
	  my $cds_point_loc = 
	    $gene_coords->map(Bio::Location::Simple->new(-seq_id=>$chrom,
							 -start=> $pos,
							 -end  => $pos));
	  #print $cds_point_loc->start, " ", $cds_point_loc->seq_id," for $chrom:$pos\n";
	  if ( ! $cds_point_loc->start ) {
	    # UTR located SNV
	    $group{'Variant_effect'} = sprintf("%s %d %s %s",
					       'sequence_variant_causing_no_change_of_translational_product',
					       $var_length,
					       $mRNA->primary_tag,
					       $name);

	    last;
	  }
	  my $n =1;
	  my %alt_CDS;

	  for my $CDS ( sort { $a->start * $a->strand <=> 
				 $b->start * $b->strand } 
			$mRNA->get_SeqFeatures('CDS') ) {
	    my ($start,$end) = ($CDS->start, $CDS->end);
	    if ( $CDS->strand < 0 ) { 
	      ($start,$end) = sort { $b <=> $a } ($start,$end);
	    }
	    $cds .= $dbh->seq($chrom,$start => $end );
	    $n++;
	  }	    
	  my $mapped_SNP = $cds_point_loc->start-1;
	  $group{'Ref_Codon'} = substr($cds,
				       3*int($mapped_SNP / 3),
				       3);
	  # iterate through the possible alleles (really there should be only 1) (Highlander!!)
	  warn("alleles = @read_alleles $chrom:$pos\n") if @read_alleles >1;
	  for my $allele ( @read_alleles ) {
	    $alt_CDS{$allele} = $cds;
	    if ( $debug ) {
	      warn("$name $chrom $pos ", $mRNA->location->to_FTstring()," $allele\n");
#	      warn("$ref_base --> $read_allele\n") if $debug;
	      warn("cds ", length($cds)," point is ", $cds_point_loc->start,"\n");
	    }
	    my $orig;
	    if ( $mRNA->strand > 0 ) { 
	      $orig = substr($alt_CDS{$allele},$mapped_SNP,1);
	    } else {
	      $orig = &reverse_complement(substr($alt_CDS{$allele},$mapped_SNP,1));
	    }
	    if ( $orig ne $ref_base ) {
	      my ($left,$right) = ($pos-1,$pos+1);
	      ($right,$left) = ($left,$right) if $mRNA->strand < 0;
	      my $base = $dbh->seq($chrom,$left => $right);
	      my $context = substr($alt_CDS{$allele},$mapped_SNP-1,3);
	      warn("orig and ref disagree ($orig:$context) ($ref_base:$base) $chrom:$pos to ($allele)\n");
	    }
	    my $allele_replace = $mRNA->strand > 0 ? $allele : &reverse_complement($allele);
	    substr($alt_CDS{$allele},$mapped_SNP,1,$allele_replace);
	    $group{'Variant_Codon'} = substr($alt_CDS{$allele},
					     3*int($mapped_SNP / 3),
					     3);
	    # figure out the codon		
	  }


	  my $pep = Bio::Seq->new(-seq=>$cds)->translate->seq;
	  $pep    =~ s/\*$//;	    
	  for my $allele ( keys %alt_CDS ) {
	    my $altpep = Bio::Seq->new(-seq => $alt_CDS{$allele})->translate->seq;		
	    if ( $debug ) {
	      open(my $t1 => ">test/$name\_$row.fa") || die $!;
	      print $t1 ">$name\n$cds\n";
	      print $t1 ">alt\n$alt_CDS{$allele}\n";
	      open(my $t2 => ">test/$name\_$row\_pep.pep") || die $!;
	      print $t2 ">$name\n$pep\n";
	      print $t2 ">alt\n$altpep\n";
	    }
	    my $effect;	    
	    $altpep =~ s/\*$//;
	    if ( $altpep ne $pep ) {		
	      if ( $altpep =~ /\*/ ) {
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
		
	  }
	  last;
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
		      $pos, $pos+1, $qual, '+',
		      '.',$group), "\n";
      last if $debug && $row > 100;

    }
  }

sub reverse_complement { Bio::Seq->new(-seq => shift)->revcom->seq }
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
