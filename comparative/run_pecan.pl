#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::TreeIO;
use IO::String;
use File::Spec;
use Getopt::Long;
my $exe = 'java -cp /usr/local/pkg/pecan/pecan.jar bp.pecan.Pecan';
my $treestring = '((n_crassa,n_tetrasperma_2508),n_discreta_8579);';

my $seqdir = 'seqs'; # working dir
GetOptions('t|tree:s' => \$treestring,
	   'wd:s'     => \$seqdir,
	   'exe:s'    => \$exe,
	  );
my $i = 0;
my $cp = $treestring;
$cp =~ s/[\(\);]//g;
my @expected = map { s/:\S+//g; $_; } split(/,/,$cp);

my $dir = shift || "pecan_alignments";
$dir = File::Spec->rel2abs($dir);
opendir(DIR, $dir) || die $!;
for my $subdir ( readdir(DIR) ) {
    next if( $subdir =~ /^\./ || ! -d "$dir/$subdir");
    chdir("$dir/$subdir");
    my $in = Bio::SeqIO->new(-format => 'fasta',
			     -file   => 'seqs.fasta');
    mkdir("$dir/$subdir/$seqdir");
    my %seqs;
    while( my $seq = $in->next_seq ) {
	my $id = $seq->display_id;
	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$dir/$subdir/$seqdir/$id");
	$out->write_seq($seq);
	$seqs{$id}++;
    }
    my $tree = Bio::TreeIO->new(-fh => IO::String->new($treestring),
				-format => 'newick')->next_tree;
    
    my @final_list;
    for my $sp ( @expected ) {
	if( ! $seqs{$sp} ){ 
	    my $node = $tree->find_node(-id => $sp);
	    $tree->remove_Node($node);
	} else {
	    push @final_list, "$seqdir/$sp";
	}
    }
    for my $node ( $tree->get_nodes ) {
	#$node->branch_length('');
    }
    my $tree_str = IO::String->new;
    my $treeout = Bio::TreeIO->new(-format => 'newick',
				   -fh     => $tree_str);
    $treeout->write_tree($tree);
    $treeout= undef;
    my $str = ${$tree_str->string_ref};
    $str =~ s/,/, /g;
    chomp($str);
    my $cmd = "$exe -E '$str' -F ".join(" ", @final_list);
    system($cmd);
    symlink("output.mfa","mavid.mfa");
    open(my $tfh => ">$dir/$subdir/job.sh") || die $!;
    print $tfh "$cmd\n";
    close($tfh);
}
