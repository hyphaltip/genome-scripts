#!/usr/bin/perl -w
use strict;
my $window = 0;
use Bio::DB::SeqFeature::Store;
use Bio::DB::Sam;
use Getopt::Long;
my $dir =  'gff3';
my $domains = 'batrachochytrium_dendrobatidis_JEL423_1.InterPro.tab';
my $bamfolder =  '.';
my $genome = ;
my $feature = 'gene:Broad';
GetOptions(
	   'd|dir:s' => \$dir,
	   'domains:s'=> \$domains,
	   'b|bamfolder:s' => \$bamfolder,
	   'g|genome:s' => \$genome,
	   'f|feature:s' => \$feature,
	   );

my $db =  Bio::DB::SeqFeature::Store->new
  (-adaptor => 'berkeleydb',
   -dir     => $dir);

my %sams;
opendir(DIR, $bamfolder) || die $!;
my @sam_names;
for my $d ( readdir(DIR)) {
  next unless $d =~ /(\S+)\.bam/;
  push @sam_names, $1;
  $sams{$1} = Bio::DB::Sam->new(-bam  => "$bamfolder/$d",
				-fasta=> $genome,
			       );
}
@sam_names = sort @sam_names;

open(my $ifh => $domains) || die $!;
my %gene_domains;
my %go;
while (<$ifh> ) {
  chomp;
  my @row = split(/\t/,$_);
  $gene_domains{$row[0]}->{join(" ","$row[3]:$row[4]",$row[5])}++;
  my $last = pop @row;
  next if $last !~ /GO/;
  $go{$row[0]}->{$_}++ for split(/,\s*/,$last);
}

my @features   = $db->get_features_by_type($feature);
open(my $fh => ">Gene_reads.tab") || die $!;
print join("\t", qw(GENE), (map { $_."_READ_COUNT" } @sam_names),
	   qw(DOMAINS GO)), "\n";

for my $f ( @features ) {
  my %count;
  for my $s ( @sam_names ) {
    my @aln = $sams{$s}->get_features_by_location(-seq_id => $f->seq_id,
						  -start  => $f->start,
						  -end    => $f->end);
    $count{$s} = scalar @aln;
  }
  print $fh join("\t",$f->name, map { $count{$_} } @sam_names,
		 join("; ", keys %{$gene_domains{$f->name} || {}}),
		 join("; ", keys %{$go{$f->name} || {}}),
		), "\n";
}
