#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $USAGE = <<EOF
orthomcl2orthonly 
REQUIRED VALUES
  -db dbfile/dir  location of fasta sequence file or directory of files
OPTIONAL OPTIONS
  --missing       number of taxa allowed to be missing in ortholog cluster (default: 0)
  --maxsp         number of gene members per species that are allowed (default: 1)
  -i mcl_file     orthomcl result file (default: all_orthomcl.out)

  -r report_file  write the ortholog report to a file, otherwise will be STDOUT
  -o  outputdir   where each ortholog sequence set should be written (default: orthomcl_orthologs)
  -p  prefix      prefix that ortholog sequence will be written with (default: orthomcl_)           

 -h/--help        display this help
EOF
    ;
    

my ($help,$dbdir,$tab,$orthomcl);
my $allowed_missing  = 0;
my $outdir = 'orthomcl_orthologs';
my $prefix = 'orthmcl_';
my $infile;
my $max_per_species = 1;
my $suffix = ".fa";

GetOptions(
	   'i|in:s'         => \$infile,
	   'db|dbdir=s'     => \$dbdir,
	   'o|out|outdir:s' => \$outdir,
	   'missing:i'      => \$allowed_missing,      
	   'maxsp:i'        => \$max_per_species,
	   'p|prefix:s'     => \$prefix,
	   'r|report:s'     => \$tab,
	   'h|help'         => \$help,
	   );
die("usage: $USAGE\n") if( $help || ! $dbdir);

mkdir($outdir) unless -d $outdir;
my $report_fh;

if( $tab ) {
    open($report_fh, ">$tab") || die("cannot open $tab: $!\n");
} else {
    $report_fh = \*STDOUT;
}
my $dbh = Bio::DB::Fasta->new($dbdir);

my $in;
open($in, $infile) || die ("cannot open $infile, provide one via -i infile: $!\n");
my @orthologs;
my %taxa;
while(<$in>) {
    my ($og, @genes) = split;
    my $n;
    if( $og =~ /(\d+)([:,.]?)$/ ) {
	$n = $1;
    } else {
	warn("unparseable group number $_ ($og)\n");
	next;
    }
    for my $gene ( @genes ) {
	my ($sp,$gname) = split(/\|/,$gene);
	$taxa{$sp} = 1;
	push @{$orthologs[$n]->{$sp}}, $gene;	
    }
}

warn(scalar @orthologs, " orthologs captured\n");
my $i = 0;
my $ntaxa = scalar keys %taxa;
warn("ntaxa = $ntaxa\n");
for my $orth ( @orthologs ) {
    my @seent = keys %{$orth};
    if( scalar @seent >= ($ntaxa - $allowed_missing) ) {	
	my @genes;
	for my $sp ( sort keys %$orth ) {
	    my $g = $orth->{$sp};
	    if( scalar @$g > $max_per_species ) {
		@genes = ();
		last;
	    } else { 
		push @genes, @$g;
	    }
	}
	if( @genes ) {	    
	    # @genes = sort @genes;
	    my $outseq = Bio::SeqIO->new(-format => 'fasta',
					 -file   => ">".File::Spec->catfile($outdir, $prefix . $i.$suffix));
	    for my $gene ( @genes) {
		my $geneseq = $dbh->get_Seq_by_acc($gene);
		if( ! defined $geneseq ) {
		    warn("cannot find $gene in db $dbdir\n");
		    next;
		}
		$outseq->write_seq($geneseq);
	    }
	    printf $report_fh "%d %s\n", $i, join("\t", @genes);
	} 	
    } else {
#	warn("saw ", scalar @seent, " taxa in orth $i\n");
    }
    $i++;
}
