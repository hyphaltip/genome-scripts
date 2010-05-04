#!/usr/bin/perl -w
use strict;
use Bio::DB::Sam;
use Getopt::Long;
use List::Util qw(sum);
use Bio::DB::SeqFeature::Store;
use Env qw(HOME);
# provide SAM file(s) and look for evidence of siRNAs 
# (common FWD and REV strand smallRNA reads in a window) 
# optionally use gene annotation as helping point 

my ($user,$pass,$host,$dbname) = ('','','localhost');

my $src = 'gene:AUGUSTUS';
my $output;
my ($min_reads) = 10;
my @bams;
my $genome;
my ($debug);
GetOptions(
	   'v|verbose!' => \$debug,
	   'u|user:s' => \$user,
	   'p|pass:s' => \$pass,
	   'host:s'   => \$host,
	   'db|dbname:s' => \$dbname,
	   'b|bam:s'   => \@bams,
	   'g|genome:s'=> \$genome,
	   's|src:s'   => \$src,
	   'm|min_reads' => \$min_reads,
	   'output:s'  => \$output,
	   );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}
unless(@bams) {
    die("need a bam file\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );

my %bam_db;
my %bam_mapped;
my @bam_names;

for my $bam ( @bams ) {
    if( $bam =~ /(\S+)\.bam/ ) {
	my $stem = $1;
	$stem =~ s/\.clp_\d+\S*//;
	$bam_db{$stem} = Bio::DB::Sam->new(-bam => $bam,
					   -expand_flags => 1,
					   -fasta => $genome);
	my $msg = `samtools flagstat $bam`;
	if( $msg =~ /(\d+)\s+mapped\s+/ ) {
	    $bam_mapped{$stem} = $1 / 1_000_000;
	} else { $bam_mapped{$stem} = 1 }
#	warn "library $bam -- ",$bam_mapped{$stem},"\n";
	push @bam_names, $stem;
    } else {
	warn("unknown BAM file $bam\n");
	next;
    }
}

if( $output && $output ne '-' ) { 
    open($output => ">$output" ) || die $!;
} else {
    $output = \*STDOUT;
}

my $iter = $dbh->get_seq_stream(-type => $src);
my (undef,$from) = split(/:/,$src);

print $output join("\t",'MRNA',
		   qw(CHROM START STOP STRAND FOUND_SIRNA),
		   ( map { $_.'.FWD_READ',
			   $_.'.FWD_RPKM',
			   $_.'.REV_READ',
			   $_.'.REV_RPKM',
			   $_.'.ALL_READ',			   
			   $_.'.ALL_RPKM',			   
		       } @bam_names ),
		   'ALL_READ_COUNT',
		   ),"\n";
while( my $gene = $iter->next_seq ) {
    my $name = $gene->name;
    my $status = 'NO';
    print join("\t", 
	       $name,
	       $gene->seq_id,
	       $gene->start,
	       $gene->end,
	       $gene->strand < 0 ? 'R' : 'F',"",
	       );
    my @row;
    my $all_total = 0;
    my $gene_len = $gene->length / 1_000;
    
    for my $bam_id ( @bam_names ) {
	my $sam = $bam_db{$bam_id};
	my %count;	
	
	for my $aln ( $sam->get_features_by_location(-seq_id => $gene->seq_id,
						     -start  => $gene->start,
						     -end    => $gene->end) ){
	    # count the number of reads, organize by strand
	    $count{$aln->strand < 0 ? 'REV' : 'FWD'}++;	
	    $count{'ALL'}++;	
	}

	my $total = (sum (values %count)) || 0;
	$all_total += $total;
	if( $total > $min_reads ) {
	    $status = 'YES';
	}
	push @row, map { ($count{$_} || 0),
			 sprintf("%.3f",
				 (($count{$_} || 0)/$gene_len)
				 / $bam_mapped{$bam_id}) } qw(FWD REV ALL);
    }
    print join("\t", $status, @row,$all_total),"\n";
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
