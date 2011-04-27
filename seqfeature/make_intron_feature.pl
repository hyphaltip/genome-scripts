#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Env qw(HOME);
my ($user,$pass,$dbname,$host);
$host ='localhost';
my $prefix;
my $debug = 0;
my $src = 'gene:NCBI_PASA';
my $output;

GetOptions(
'v|verbose!' => \$debug,
'u|user:s' => \$user,
'p|pass:s' => \$pass,
'host:s' => \$host,
'db|dbname:s' => \$dbname,

's|src:s' => \$src,
'o|output:s' => \$output,
);

unless( defined $dbname ) {
    die("no dbname provided\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn => $dsn,
                                          -user => $user,
                                          -password => $pass,
                                          );
my $ofh;
if( $output && $output ne '-' ) {
    open($ofh => ">$output" ) || die $!;
} else {
    $ofh = \*STDOUT;
}

print $ofh "##gff-version 3\n";
my $iter = $dbh->get_seq_stream(-type => $src);
my (undef,$from1) = split(/:/,$src);
my $count = 0;
my $gene_count = 0;

while( my $gene = $iter->next_seq ) {
    $gene_count++;
    my $gene_name = $gene->name;
    my @mRNA = $gene->get_SeqFeatures('mRNA');
    if( ! @mRNA ) {
	warn("no mRNA for $gene_name\n");
    }
    for my $mRNA ( @mRNA ) {	# 1st mRNA for now
	my $last_exon;
	my $i = 1;
	my $mRNA_name = $mRNA->name || 'mRNA-'.$mRNA->load_id;
	for my $exon ( sort { $a->start * $a->strand <=> $b->start * $b->strand }
		       $mRNA->get_SeqFeatures('exon') ) {
	    if( $last_exon ) {
		my ($start,$end) = ( $last_exon->end+1,$exon->start - 1);
		if( $exon->strand < 0 ) {
		    ($start,$end) = ( $exon->end + 1, $last_exon->start-1);
		}

		print $ofh join("\t",
				$gene->seq_id,
				$gene->source,
				'intron',
				$start,$end, '.',
				$exon->strand > 0 ? '+' : '-',
				'.',
				sprintf('ID=%s.i%d;Gene=%s',
					$mRNA_name,
					$i++,
					$gene_name)),"\n";
	    }
	    $last_exon = $exon;
	}
	last;
    }
    last if $debug && $count++ > 10;
}

warn("gene count $gene_count\n");
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

