#!/usr/bin/perl -w
use strict;
use Bio::TreeIO;
use Statistics::Descriptive;
use Getopt::Long;

my $master_tree;
my $dir = '.';
my $extension = '.tre';
GetOptions('t|tree:s' => \$master_tree,
	   'd|dir:s'  => \$dir,
	   'e|ext:s'  => \$extension);


opendir(DIR,$dir)|| die $!;

my $in = Bio::TreeIO->new(-format => 'newick',
			  -file   => $master_tree);
my %branches;
my $t = $in->next_tree;
if( ! $t ) { die "cannot find a tree in file $master_tree\n"; }

for my $node ( $t->get_nodes ) {
    next unless defined $node->ancestor;
    my $id = $node->id;
    if( ! $id ) {
	$id = join(",", sort map { $_->id } grep { $_->is_Leaf } $node->get_all_Descendents);
    }
    my $ancestor = $node->ancestor;
    my $anc_id;
    if( ! $ancestor ) {
	$anc_id = 'root';
    } else {
	$anc_id = join(",", sort map { $_->id } grep { $_->is_Leaf } $ancestor->get_all_Descendents);
    }
    $branches{"$id,$anc_id"} = Statistics::Descriptive::Full->new();
}

for my $file ( readdir(DIR) ) {

    next unless $file =~ /\Q$extension\E$/;
    my $in = Bio::TreeIO->new(-format => 'newick',
			      -file   => "$dir/$file");
    if( my $tree = $in->next_tree ) {
	for my $node ( $tree->get_nodes ) {
	    next unless defined $node->ancestor;
	    my $id = $node->id;
	    if( ! $id ) {
		$id = join(",", sort map { $_->id } grep { $_->is_Leaf } $node->get_all_Descendents);
	    }
	    my $ancestor = $node->ancestor;
	    my $anc_id;
	    if( ! $ancestor ) {
		$anc_id = 'root';
	    } else {
		$anc_id = join(",", sort map { $_->id } grep { $_->is_Leaf } $ancestor->get_all_Descendents);
	    }
	    if( ! exists $branches{"$id,$anc_id"} ) {
		warn("cannot find node for  $id -> $anc_id branch in $file\n");
		next;
	    }
	    $branches{"$id,$anc_id"}->add_data($node->branch_length);
#	    warn("$id -> $anc_id = ",$node->branch_length || 0,"\n");
	}
    }
    last;
}


for my $node ( $t->get_nodes ) {
    next unless defined $node->ancestor;
    my $id = $node->id;
    if( ! $id ) {
	$id = join(",", sort map { $_->id } grep { $_->is_Leaf } $node->get_all_Descendents);
    }
    my $ancestor = $node->ancestor;
    my $anc_id;
    if( ! $ancestor ) {
	$anc_id = 'root';
    } else {
	$anc_id = join(",", sort map { $_->id } grep { $_->is_Leaf } $ancestor->get_all_Descendents);
    }
    $node->branch_length($branches{"$id,$anc_id"}->mean);
}
my $out = Bio::TreeIO->new(-format => 'newick',
			   );
$out->write_tree($t);
