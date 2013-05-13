#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Bio::DB::Fasta;

my ($input,$db,$odir);
my $ext = 'pep';
my $min_size = 1;
my $orthomcl = 1;
GetOptions(
	   'i|input:s'  => \$input,
	   'o|output:s' => \$odir,
           'orthomcl!'  => \$orthomcl,
	   'd|db:s'     => \$db,
           'ext:s'      => \$ext,
	   'min:i'      => \$min_size, 
	   );

$input ||= shift @ARGV;
$odir ||= $input.'-seqs';

my $dbh = Bio::DB::Fasta->new($db);
open(my $fh => $input) || die "cannot open $input: $!\n";
mkdir($odir) unless -d $odir;
my $groupct = 1;
while(<$fh>) {   
    my ($group,@orthologs);
    if( $orthomcl ) {
	($group,@orthologs) = split;    
	$group =~ s/://;
    } else {
	@orthologs = split;
	$group = sprintf("FAM_%05d",$groupct++);
    }
    next if( @orthologs < $min_size );
    my %count;
    my %domains; 
    my $seqio =Bio::SeqIO->new(-format => 'fasta', 
			       -file => ">$odir/$group.$ext.fa");
    for my $orth (@orthologs ) {
     my $seq = $dbh->get_Seq_by_acc($orth);
     if( ! $seq ) { warn("cannot lookup $orth in the DB!"); }
     else { 
       $seqio->write_seq($seq);
     }
     #$seqio->write_seq( map { $dbh->get_Seq_by_acc($_) } @orthologs);
    }
}
