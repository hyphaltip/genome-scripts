#!/usr/bin/perl -w
use strict;

use Env qw(HOME);
use Getopt::Long;
use Bio::DB::SeqFeature;
my $debug = 0;

my $method  = 'predicted';
my $host   = 'localhost';
my ($user,$pass);
GetOptions('v|verbose|debug!' => \$debug,
	   'm|method:s'       => \$method,
	   'u|user:s'         => \$user,
	   'p|pass|passwd:s'  => \$pass,
	   );

my $dir = shift;
my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'berkeleydb',
					 -dir     => $dir);

print "#", join("\t", qw(gene chrom strand start stop source method)),
    "\n";

my %CHROMS;
{
    my $i = 0;
    for my $chrom ( sort { $a->id cmp $b->id } 
		    $db->get_features_by_type('scaffold:chromosome') ) {
	$CHROMS{$chrom->seq_id} = [$i++,$chrom->length];
    }
}

for my $chrom ( sort { $CHROMS{$a}->[0] <=> $CHROMS{$b}->[0] }
		keys %CHROMS ) {
    my (@segments) = $db->segment(-name => $chrom,
				  -type => "scaffold:chromosome");
    warn( "segments are @segments\n") if @segments > 1;
    my $segment = shift @segments;
    for my $gene ( sort { $a->start * $b->strand <=> 
			      $b->start * $b->strand } 
		   $segment->features('tRNA_gene:tRNAscan-SE') ) {
	print join("\t", $gene->id, $gene->seq_id,
		   $gene->strand > 0 ? "+1" : "-1", 
		   $gene->start, $gene->end, 
		   $gene->source_tag, $method),"\n";
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
