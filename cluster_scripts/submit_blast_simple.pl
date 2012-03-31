#!/usr/bin/perl -w
use strict;
use File::Spec;
use Bio::SeqIO;
use Getopt::Long;

my $cpus = 4;
my $home = $ENV{'HOME'};
my $jobdir = File::Spec->catfile($home,'jobs');
my $input  = File::Spec->catfile($home,'input');
my $outdir = File::Spec->catfile($home,'output');
my $tmpdir = File::Spec->catfile($home, 'tmp');
my $qsize = 5000;
my ($query,$dbin,$params,$exe);
$exe = 'blastall'; # default
my $prefix = 'blast';
$params = '-p blastp -e 1e-3 -m 8';
my $uid;
GetOptions(
	   'i|q|query:s' => \$query,
	   'd|db:s'      => \$dbin,
	   'p|params:s'  => \$params,
	   's|size:i'    => \$qsize,
	   'e|exe:s'     => \$exe,
           'id|uid:s' => \$uid,
           'c|cpus:i'       => \$cpus,
	   'indir|input:s'  => \$input,
	   'o|out|outdir:s' => \$outdir,
	   'j|job|jobdir:s' => \$jobdir,
	   't|tmp|tmpdir:s' => \$tmpdir,
	   );
die( "need a dbname\n" ) unless (defined $dbin );
my @dblst;
for my $dbpart ( split(/\s+/,$dbin) ) {
    my $rel = File::Spec->rel2abs($dbpart);
    my (undef,undef,$dbname) = File::Spec->splitpath($rel);
    push @dblst, $rel;# File::Spec->catfile($scratchdir,$prefix,$dbname);    
} 
my $db = join(" ", @dblst);
warn("db is $db\n");
die( "need an inputfile\n" ) unless defined $query && -e $query;
die( "need some params\n" ) unless defined $params;


for my $dir( $jobdir,$input,$outdir,$tmpdir ) {
    mkdir $dir unless -d $dir;
}

$query = File::Spec->rel2abs($query);
my (undef,undef,$query_file) = File::Spec->splitpath($query);
$uid = $query_file unless defined $uid;
my $in = Bio::SeqIO->new(-file   => $query,
			 -format => "fasta");

my ($ct,$setct) = (0,1);
my $infile = File::Spec->catfile($input,"$uid.$setct.fa");
my $input_seqfh = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$infile");
				  
my @infiles = ($infile);
my @infilenames = ("$uid.$setct.fa");

while( my $s = $in->next_seq ) {
	next if $s->seq =~ /^[Xx]+$/;
    $input_seqfh->write_seq($s);
    if( ++$ct >= $qsize ) {
	$input_seqfh->close();
	$setct++;
	$infile = File::Spec->catfile($input,"$uid.$setct.fa");
	$input_seqfh = Bio::SeqIO->new(-format => 'fasta',
				       -file   => ">$infile");
	push @infiles, $infile;
	push @infilenames, "$uid.$setct.fa";
	$ct = 0;
    }
}

my $jobfile = ">$jobdir/$prefix\_$uid.sh";
open my $fh => ">$jobfile" || die $!;
chmod(0700,$jobfile);
print $fh ( "#!/bin/bash\n",
	    '#PBS'." -N $uid\n",
	    "#PBS -l nodes=1:ppn=$cpus\n",
	    "source ~/.bash_profile\n",
    );
print $fh "$exe $params -a $cpus -i $input/$uid.\$PBSARRAY_ID.fa -d \"$db\" -o $outdir/$uid.\$PBSARRAY_ID.tab \n";


print "qsub -t 1-$setct%50 -j oe -o $outdir/$uid.BLAST.out $jobfile\n";
