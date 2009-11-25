#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Env qw(HOME);
my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $other = 'gene:JGI';
my $src = 'gene:GLEAN_PASA';
my $output;

GetOptions(
	   'v|verbose!' => \$debug,
	   'u|user:s' => \$user,
	   'p|pass:s' => \$pass,
	   'host:s'   => \$host,
	   'db|dbname:s' => \$dbname,
	   
	   's|src:s'   => \$src,
	   'o|other:s' => \$other,
	   
	   'output:s'  => \$output,
	   );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

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

my $iter = $dbh->get_seq_stream(-type => $src);
my (undef,$from1) = split(/:/,$src);
my (undef,$from2) = split(/:/,$other);
print $output join("\t", $from1,'MRNA_FROM',qw(CHROM START_FROM STOP_FROM STRAND_FROM), 
		   $from2, 'MRNA_TO', qw(START_TO STOP_TO STRAND_TO SINGLE_MATCH)), "\n";
while( my $gene = $iter->next_seq ) {
    my $name = $gene->name;
    my ($mRNA) = $gene->get_SeqFeatures('mRNA'); # 1st mRNA for now
    if( my @overlaps = $gene->segment->features(-type => $other) ) {
	for my $overlap ( @overlaps ) {
	    my ($mRNA_overlap) = $overlap->get_SeqFeatures('mRNA');
	    print $output join("\t", $name, $mRNA->name, $gene->seq_id,
			       $gene->start,$gene->stop,$gene->strand,
			       $overlap->name, $mRNA_overlap->id,
			       $overlap->id, $overlap->start,$overlap->end,
			       $overlap->strand, @overlaps == 1 ? 'yes' : 'no'), "\n";
	}
    } else {
	print $output join("\t", $name, $gene->seq_id,
			       $gene->start,$gene->stop,$gene->strand,
			   '','','','','NO_MATCH'),"\n";
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
