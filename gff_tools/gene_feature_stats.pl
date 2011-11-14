#!/usr/bin/perl -w
# $Id$

=head1 NAME

gene_feature_stats - summary of gene features in a GFF file

=head1 USAGE

gene_feature_stats gffdir

=head1 DESCRIPTION 

CMDLINE options
 -genome [size] genome size to give fraction of genome

=cut

use strict;
use Getopt::Long;
use Statistics::Descriptive;
use Bio::DB::SeqFeature;
use Env qw(HOME);

my ($user,$password,$dbname);
my $host = 'localhost';
my $dir = shift;
GetOptions("db|dbname:s" => \$dbname,
	   "u|user:s"    => \$user,
	   "p|pass|password:s" => \$password);

($user,$password) = &read_cnf unless defined $user;

#my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$dbname,$host);

my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'berkeleydb',
					  -dir     => $dir);
#my  $db = Bio::DB::SeqFeature::Store->new(-adaptor  => 'DBI::mysql',
#					  -user     => $user,
#					  -password => $password,
#					  -dsn      => $dsn);
my $genome = 0;
my %CHROMS;
{
    my $i = 0;
    for my $chrom ( sort { $a->id cmp $b->id } 
		    $db->get_features_by_type('scaffold') ) {
	$CHROMS{$chrom->seq_id} = [$i++,$chrom->length];
	$genome += $chrom->length;
    }
}

my $total;
my $count;
my %stats;
my @ftypes = qw(gene mrna cds exon five_prime_utr three_prime_utr intron);
for my $t ( @ftypes ) {
    $stats{$t} = Statistics::Descriptive::Full->new;
}

for my $chrom ( sort { $CHROMS{$a}->[0] <=> $CHROMS{$b}->[0] }
		keys %CHROMS ) {
    my $segment = $db->segment($chrom);
    my $lastgene;
    
    for my $gene ( $segment->features('gene') ) {
	$stats{'gene'}->add_data($gene->length);
	for my $mRNA( $gene->get_SeqFeatures ) {
	    for my $t ( qw(exon cds) ) {
		for my $e ( $mRNA->get_SeqFeatures($t) ) {
		    $stats{$t}->add_data($e->length);
		}
	    }
	    my $mRNAlen = 0;
	    for my $e ( $mRNA->get_SeqFeatures('exon') ) {
		$mRNAlen += $e->length;
	    }
	    $stats{'mrna'}->add_data($mRNAlen);
	    for my $t ( qw(five_prime_utr three_prime_utr) ) {
		my $utr_len = 0;
		for my $utr ( $mRNA->get_SeqFeatures($t) ) {
		    $utr_len += $utr->length;
		}
		$stats{$t}->add_data($utr_len) if $utr_len;
	    }
	    my $lastexon;
	    for my $exon ( sort { $a->start * $a->strand <=> 
				      $b->start * $b->strand  }
			   $mRNA->get_SeqFeatures('exon') ) {
		if( $lastexon ) {
		    my $intronlen;
		    if( $exon->strand < 0 ) {
			$intronlen = $lastexon->start - ($exon->end+1);
		    } else {
			$intronlen = $exon->start - ($lastexon->end+1);
		    }
		    $stats{'intron'}->add_data($intronlen);
		}
		$lastexon = $exon;
	    }
	}    
    }
}

for my $t (  @ftypes ) {
    printf "%s Max = %d Mean = %.1f Median = %1.f Total = %d N = %d %.2f Mb Total %% %.2f\n",
    $t, $stats{$t}->max, $stats{$t}->mean, $stats{$t}->median,
    $stats{$t}->sum,$stats{$t}->count,
    $stats{$t}->sum / 1000000, 
    100 * $stats{$t}->sum / $genome;
}

print "Genome size is $genome\n";

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
