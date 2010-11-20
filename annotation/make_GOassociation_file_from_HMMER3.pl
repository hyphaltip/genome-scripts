#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use POSIX;
use GO::AppHandle;
use Env;
    
my $srcPrefix = 'MINE';
my $centerPrefix = 'MYLAB';
my $minE = 1;
my $taxon = '0';
# expect HMMER3 --domtabl TABLE results and PFAM2GO file
my $dbname = 'GO';
GetOptions('c|p|prefix:s'  => \$centerPrefix,
	   't|taxon:s' => \$taxon,
	   'm|minE:f'  => \$minE,
	   's|src:s'   => \$srcPrefix,	   
	   'db|dbname:s'    => \$dbname,
	   );

my ($user,$pass,$dsn) = read_cnf();
$user ||= $USER;
$dsn ||= 'localhost';
#warn("dbname is $dbname\n");
my $apph = GO::AppHandle->connect(-dbname =>$dbname, 
				  -dbhost =>$dsn,
				  -dbuser => $user,
				  -dbauth =>$pass,);


my ($hmmertable, $pfam2go_f) = @ARGV;
die("need hmmertable and pfam2go filenames as input\n") 
    unless defined $hmmertable && defined $pfam2go_f;

my %pfam2go;
if( $pfam2go_f =~ /\.gz$/) {
        open(PFAM, "zcat $pfam2go_f |") || die "$pfam2go_f: $!";
} else {
    open(PFAM, $pfam2go_f) || die "$pfam2go_f: $!";
}
while(<PFAM>) {
    next if /^\s*\!/;
    chomp;
    next if ( /^\!/);
    my ($pfamstr, $goset) = split(/\s+>\s+/,$_);
    my ($pfamid, $name) = split(/\s+/,$pfamstr);
    my ($godesc,$goid) = split(/\s+;\s+/,$goset);
    push @{$pfam2go{$name}}, [$goid,$godesc];
}
close(PFAM);
my %genome2pfam;
my $date = strftime("%d-%b-%Y", localtime);
my %seen;
open(HMMERTABLE, $hmmertable) || die "$hmmertable: $!";
while(<HMMERTABLE>){
    chomp;
    next if /^\#/ || /^\s+$/;
    my ($domain,$domacc,$tlen,$qname,$qacc,$qlen, $fullevalue,$fullscore,$fullbias,
	$n,$ntotal,$cvalue,$ivalue,$domscore,$dombias,$hmmfrom,$hmmto,$alifrom,$alito, $envfrom,$envto) = split(/\s+/,$_,23);
    next if $fullevalue > $minE;
    if( $pfam2go{$domain} ) {
	for my $go ( @{$pfam2go{$domain}} ) {
	    my ($goid,$godesc) = @$go;
	    my $term = $apph->get_term({acc => $goid});
	    my (undef,$code) = split(/_/,$term->type);
	    $code = uc substr($code,0,1);
	    next if $seen{$qname."_".$goid}++;
	    print join("\t", 
		       $srcPrefix,
		       $qname, $qname, '',
		       $goid, 'HMMScan', 'ISS', '', $code,
		       '','', 'protein', "taxon:$taxon",
		       $date,
		       $centerPrefix),"\n";
	}
    }
}
close(HMMERTABLE);


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
    
