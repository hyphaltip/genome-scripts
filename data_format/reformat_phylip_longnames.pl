#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::AlignIO;

my $new_id_len = 20;
my $ext = 'nex';
GetOptions(
	   'l|len:s' => \$new_id_len,
	   'ext:s'   => \$ext,
	   );

my $dir = shift || ".";

opendir(DIR, $dir) || die $!;
for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.\Q$ext\E$/;
    my $stem = $1;
    my $in = Bio::AlignIO->new(-format => 'nexus',
			       -file   => "$dir/$file");
    my $out = Bio::AlignIO->new(-format => 'phylip',
				-file   => ">$dir/$stem.phy",
				-interleaved => 0,
				-idlength=> $new_id_len);
    while( my $aln = $in->next_aln ) {
	$out->write_aln($aln);	
    }
}
