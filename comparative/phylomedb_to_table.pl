#!/usr/bin/perl -w
use strict;

# query is the file which has the names of genes you want to lookup. Typically the ORF names like YAL001C but CHS1 will also work if alias is there
# query_sp is the name of the species (5letter abbreviation like YEAST or SCHPO) that you want to report on
# gene_table is the name of the phylomedb all_gene_names.txt gene table which you can get here ftp://phylomedb.org/phylomedb/all_gene_names.txt.gz
# phylome_file is the name of the phylome you want to work with, Phylome3 is the yeast P60 dataset (60 fungi with Yeast as focus)
#  available here ftp://phylomedb.org/phylomedb/phylomes/phylome_0003/phylome_0003_orthologs.tar.gz - after this is untarred you want the phylome_3.txt file

# assuming input files are uncompressed, get fancy here later if you want
# to use dynamic uncompressing
my ($query, $query_sp, $gene_table, $phylome_file) = @ARGV;

my @query_sp;
open(my $qfh => $query_sp) || die $!;
while(<$qfh>) {
 chomp;
 push @query_sp, $_;
}
my (%genes,%order);
open(my $fh => $query) || die $!;
while (<$fh>) {
  my ($gene) = split;
  $order{$gene} = keys %order;
  $genes{$gene} = 1;
}
my %allgenes;
my %id2gene_targets;
open($fh => $gene_table) || die $!;
while (<$fh>) {
  my ($lookup,$gene) = split;
  if ( $genes{$gene} ) {
    $id2gene_targets{$lookup} = $gene;
   # simply let last writer win for now
  }
  $allgenes{$lookup} = $gene if ! exists $allgenes{$lookup};
}
warn("saw in ", scalar keys %allgenes, " genes\n");
open($fh => $phylome_file) || die $!;
my %orthologs;
while (<$fh>) {
  next if /^\#/;
  my ($name,$orthoid,$style,$CS,$trees,@coorthologs) = split;
  my ($pref,$sp) = split('_',$name);
  $pref =~ s/Phy//;
  #warn("$name -> $pref -- $sp\n");
  my ($opref,$osp) = split('_',$orthoid);
  if ( my $fullname = $id2gene_targets{$pref}) {
    push @{$orthologs{$fullname}->{$osp}}, $opref;
#    if ( @coorthologs ) {
#      for ( @coorthologs ) {
#	$orthologs{$_}->{$orthoid}++;
#      }  
  }
}
print join("\t", qw(GENE),@query_sp),"\n";
for my $orth ( sort {$order{$a} <=> $order{$b} } keys %orthologs ) {
  print join("\t", $orth, map { exists $orthologs{$orth}->{$_} ? scalar @{$orthologs{$orth}->{$_}} : 0 } @query_sp),"\n";
}


