#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;

# Algorithm
# -- get coordinates of alignment blocks,
# -- dump feature GFF in new coordinate space 

use Bio::DB::SeqFeature::Store;
use Env qw(USER HOME);
use Getopt::Long;

my $mapfile;

my ($user,$pass,$dbname,$host);
$host ='localhost';
my $outfile;
my $debug = 0;
$dbname = 'gb_neurospora_crassa_reannotation';
my $Gff = 0;
GetOptions(
	   'v|verbose!' => \$debug,
	   'u|user:s' => \$user,
	   'p|pass:s' => \$pass,
	   'host:s'   => \$host,
	   'db|dbname:s' => \$dbname,
	   'g|gff!'   => \$Gff, # GFF or BED
	   'o|out|output:s' => \$outfile,
	   'm|mapfile:s' => \$mapfile,
	   );

unless(  defined $dbname ) {
    die("no dbname provided\n");
}

unless( defined $mapfile ) {
    die("no mapfile provided\n");
}

($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);
my $dbh = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
                                          -dsn     => $dsn,
                                          -user    => $user,
                                          -password => $pass,
                                          );

open(my $mapfh => $mapfile ) || die "$mapfile: $!";
my $ofh;
if( $outfile ) {
    open($ofh => ">$outfile") || die $!;
} else {
    $ofh = \*STDOUT;
}
my $ref_genome_index = 0;
my @aln_lookups;
while(<$mapfh>) {
    my ($aln_id, @line) = split;
    my ($chrom,$start,$end) = map { $line[$ref_genome_index + $_] } 0..2; 
    
    next if $chrom eq 'NA';
    my $segment = $dbh->segment($chrom, $start => $end);
    
    for my $gene ( $segment->features('gene:NCBI_PASA') ) {	
	&write($ofh,$aln_id,$start,$gene);
	for my $mRNA ( $gene->get_SeqFeatures ) {
	    &write($ofh,$aln_id,$start,$mRNA);
	    my $i = 0;
	    for my $f (  sort { $a->start * $a->strand <=> 
				    $b->start * $b->strand } 
			 $mRNA->get_SeqFeatures('five_prime_utr') ) {
		$f->name(sprintf("%s.5utr%d",$mRNA->name,$i++));
		&write($ofh,$aln_id,$start,$f);
	    }
	    for my $f (  sort { $a->start * $a->strand <=> 
				    $b->start * $b->strand } 
			 $mRNA->get_SeqFeatures('CDS') ) {
		$f->name(sprintf("%s.cds%d",$mRNA->name,$i++));
		&write($ofh,$aln_id,$start,$f);
	    }	    
	    for my $f (  sort { $a->start * $a->strand <=> 
				    $b->start * $b->strand } 
			 $mRNA->get_SeqFeatures('three_prime_utr') ) {
		$f->name(sprintf("%s.3utr%d",$mRNA->name,$i++));
		&write($ofh,$aln_id,$start,$f);
	    }
	}
	last if $debug;
    }
    last if $debug;
}

sub write {
    if( $Gff ) {
	&write_gff(@_);
    } else {
	&write_bed(@_);
    }
}
sub write_bed {
    my ($fh,$chrom,$offset,$feature) = @_;
    my $name = $feature->name || $feature->id || $feature->display_id;

    return if ( $feature->start < $offset );
    print $fh join("\t", 
		   $chrom,
		   $feature->start - $offset, # bed format
		   $feature->end - $offset,
		   $name),"\n";
}

sub write_gff {
    my ($fh,$feature,$offset) = @_;
    my @group = ( sprintf("ID=%s",$feature->id),
		  sprintf("Name=%s",$feature->name || $feature->id));

    for my $tag ( qw(ID Parent Note Alias Name ) ) {
	if( $feature->has_tag($tag) ) {
	    my $v = join(",",$feature->get_tag_values($tag));
	    if( $v =~ /\s+/ ) {
		$v = "\"$v\"";
	    }
	    push @group, sprintf("%s=%s", $tag, $v);
	}
    }
    print $fh join("\t", 
		   $feature->seq_id,
		   $feature->source,
		   $feature->type,
		   $feature->start - $offset + 1,
		   $feature->end - $offset +1,
		   $feature->score || '.',
		   $feature->strand > 0 ? '+' : '-',
		   join(";", @group)),"\n";
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

