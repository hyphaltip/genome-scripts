#!/usr/bin/perl -w
use strict;
use Env qw(HOME USER);
use lib "$HOME/src/bioperl/core";
use Bio::TreeIO;
use Bio::TreeIO::TreeEventDBBuilder;
use Bio::DB::Tree::Store;

my $dir = shift;
my $burnin = 100; # how many trees to skip from the beginning of each file
my $dbh = Bio::DB::Tree::Store->new(-adaptor=>'DBI::SQLite',
				    -create => 1,
                                    -dsn    => 'dbname=mbtrees.idx');

my $dbh_keep = Bio::DB::Tree::Store->new(-adaptor=>'DBI::SQLite',
					 -create => 1,
					 -dsn    => 'dbname=mbtrees_keep.idx');

my $handler = Bio::TreeIO::TreeEventDBBuilder->new(-store => $dbh);
opendir(DIR, $dir) || die "$dir: $!";

for my $file ( readdir(DIR) ) {
  next unless $file =~ /\.t$/;
  warn "file is $file\n";
  my $treeio = Bio::TreeIO->new(-format  => 'nexus',
				-file    => $file,
				-handler => $handler);
  warn($treeio->_eventHandler,"\n");

#  my $treeio = Bio::TreeIO->new(-format  => 'nexus',
#				-file    => $file);
  while ( my $t = $treeio->next_tree ) {
    print "saw a tree\n";
    last;
  }
}
