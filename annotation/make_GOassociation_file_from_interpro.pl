#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use POSIX; # for strftime

use GO::Parser;
use Bio::DB::Taxonomy;
use Env;

my %uncompress = ('bz2' => 'bzcat',
		  'gz'  => 'zcat',
		  );

my $taxdb = Bio::DB::Taxonomy->new(-source => 'entrez');    
my $srcPrefix = 'UCR';
my $centerPrefix = 'StajichLab';
my $taxon;
my $species;
my $interprotab;
my $outputfile;
my $use_stdout = 0;
my $strip = 1;
my $debug = 0;
my $gofile;
GetOptions('ctr:s'       => \$centerPrefix,
	   'src:s'       => \$srcPrefix,
	   't|taxonid:s' => \$taxon,
	   's|species:s' => \$species,
	   'i|input:s'   => \$interprotab,
	   'o|output:s'  => \$outputfile,
	   'g|go:s'      => \$gofile,
	   'stdout!'     => \$use_stdout,
	   'strip!'      => \$strip,
	   'v|debug'     => \$debug,
	   );
if( $species ) {
    $taxon = $taxdb->get_taxonid($species);
}
if( ! $taxon ) {
    warn("need a taxonid, either provide it with --taxonid or give a --species='genus species' [strain] (add strain if specified in genbank)");
    exit;
}
unless($interprotab ) {
    $interprotab = shift @ARGV if @ARGV;
}

unless($interprotab && -f $interprotab ) {
    warn("must specify an interpro tab file either on the cmdline or with the -i option\n");
    exit;
}

unless( $gofile ) {
    warn("Need a GO file (gene_ontology.obo)\n");
    exit;
}
my $parser = GO::Parser->new({handler => 'obj',use_cache=>1});
#$parser->handler->file("output.xml");
$parser->parse($gofile);
my $go_graph = $parser->handler->graph;
warn("done parsing\n");


my ($stem,$fh);
if( $interprotab =~ /(\S+)\.(gz|bz2)$/) {
    $stem = ($1);
    open($fh, "$uncompress{$2} $interprotab |") || die "$uncompress{$2} $interprotab: $!";
} else {
    open($fh, "< $interprotab") || die "$interprotab: $!";
    $stem = $interprotab;
}

my $ofh;
if($use_stdout ) {
    $ofh = \*STDOUT;
} else {
    if( $strip ) {
	$stem =~ s/\.tab$//;
	$stem =~ s/\.InterPro$//i;
    }
    open($ofh,">$stem.InterPro.association");
}

my %seen;
my $date = strftime("%d-%b-%Y", localtime);
while(<$fh>) {
    # columns from IPRscan are
    chomp;
    my ($gene,$crc64,$length,$analyis_method,
	$hit_dbid,$hitdesc,$qstart,$qend,$evalue,$hit_status,
	$date_run,$iprdomain,$iprdesc,
	$go_info) = split(/\t/,$_);
    $gene =~ s/\.T\d+$//;
    if($go_info) {
	for my $goid ( split(/\|/,$go_info) ) {
	    unless( $goid =~ /(GO:\d+)/ ) {
		warn("No GOID in $goid\n");
		next;
	    }
	    next if $seen{$gene}->{$goid}++;
	    my $term = $go_graph->get_term($goid);
	    my $code;
	    if( ! $term ) {
	     $code = 'U';
            } else {
	     $code = $term->get_code_from_namespace;
	    }
	    print $ofh join("\t", 
			    $srcPrefix,
			    $gene, $gene, '',
			    $goid, $analyis_method, 'ISS', '', $code,
			    '','', 'protein', "taxon:$taxon",
			    $date,
			    $centerPrefix),"\n";
	}

	last if $debug;
    }
}
close($fh);

unless( $use_stdout ) {
    close($ofh);
}


=head1 annotation association format

    Column  Cardinality   Contents          
    ------  -----------   -------------------------------------------------------------
        0       1         Database abbreviation for the source of annotation (eg SGD)
        1       1         Database identifier of the annotated entity
        2       1         Standard name of the annotated entity
        3       0,1       NOT (if a gene is specifically NOT annotated to the term)
        4       1         GOID of the annotation     
        5       1,n       Reference(s) for the annotation 
        6       1         Evidence code for the annotation
        7       0,n       With or From (a bit mysterious)
        8       1         Aspect of the Annotation (C, F, P)
        9       0,1       Name of the product being annotated
       10       0,n       Alias(es) of the annotated product
       11       1         type of annotated entity (one of gene, transcript, protein)
       12       1,2       taxonomic id of the organism encoding and/or using the product
       13       1         Date of annotation YYYYMMDD
       14       1         Assigned_by : The database which made the annotation

=cut
    
