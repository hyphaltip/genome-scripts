#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use Bio::AlignIO;

my $debug = 0;
my $outformat = 'stockholm';
my ($alndir,$odir);
my $stem = 'mavid.mfa';
my $informat = 'fasta';
GetOptions('a|aln|dir:s'   => \$alndir,
	   'o|out:s'       => \$odir,
	   's|stem:s'      => \$stem,
	   'if|informat:s' => \$informat,
	   'of|outformat:s'=> \$outformat,
	   );


die" must provide a valid aln dir with -a or --aln or --dir\n" unless( defined $alndir);
die" must provide a valid out dir with -o or --out\n" unless( defined $odir);

mkdir($odir) unless -d $odir;
opendir(DIR, $alndir) || die $!;

my $mapfile = File::Spec->catfile($alndir,'map');
open(my $fh => $mapfile) || die $!;
my @map;
while(<$fh>) {
    chomp;
    my ($id,@row) = split(/\t/,$_);
    $map[$id] = \@row;
}

for my $dir ( readdir(DIR) ) {
    next if $dir =~ /^\./;
    my $target = File::Spec->catfile($alndir,$dir,$stem);
    next unless -f $target && ! -z $target;
    my $in = Bio::AlignIO->new(-format => $informat,
			       -file   => $target);
    my $ofile = File::Spec->catfile($odir,"$dir.stk");
    my $out = Bio::AlignIO->new(-format => $outformat,
				-file   => ">$ofile");
    if( my $aln = $in->next_aln ) {
	$aln->set_displayname_flat(1);
	$out->write_aln($aln);
    }
    last if $debug;
 }

