#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Bio::DB::Fasta;

my ($input,$db,$odir);

GetOptions(
	   'i|input:s'  => \$input,
	   'o|output:s' => \$odir,
	   'd|db:s'     => \$db,
	   );

$input ||= shift @ARGV;
$odir ||= $input.'-seqs';

my $dbh = Bio::DB::Fasta->new($db);
open(my $fh => $input) || die "cannot open $input: $!\n";
mkdir($odir) unless -d $odir;
while(<$fh>) {   
    my ($group,@orthologs) = split;    
    $group =~ s/://;
    my %count;
    my %domains; 
    my $seqio =Bio::SeqIO->new(-format => 'fasta', 
			       -file => ">$odir/$group.pep.fa");

    $seqio->write_seq( map { $dbh->get_Seq_by_acc($_) } @orthologs);
}
