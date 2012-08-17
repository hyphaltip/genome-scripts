#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqIO;
my $ext = '\.domtbl\.tab';
my $cutoff = 1e-4;
my $pfam_seqs_out = 'pfam_seq_out';
my $dbdir;
GetOptions('c|cutoff:s' => \$cutoff,
	   'e|ext:s'    => \$ext,
	   'db:s' => \$dbdir,
	   'o|out:s'=> \$pfam_seqs_out);

my $dir = shift || 'pfam';

my $dbh;
if( $dbdir && (-f $dbdir || -d $dbdir) ) {
    $dbh = Bio::DB::Fasta->new($dbdir);
    mkdir($pfam_seqs_out);
}
opendir(DIR,$dir) || die "$dir: $!";
my (%domains,%spnames,%domain_seqs);
for my $file ( readdir(DIR) ) {
    if( $file =~ /(\S+)$ext/ ) {
	my $sp = $1;
	$spnames{$sp}++;
#	print "$1 $file\n";
	open(my $fh => "$dir/$file" ) || die $!;
	while(<$fh>) {
	    next if /^\#/;
	    chomp;
	    my ($domain,$acc,$tlen,$gene,$trans,$qlen,$evalue_full,$score_full,$bias_full,
		$n, $tot,$c_evalue,$i_evalue, $domain_score,$domain_bias,
		$hmm_start,$hmm_end, $ali_start, $ali_end, 
		$env_start, $env_end,
		$accuracy, 
		$description) = split(/\s+/,$_,24);
	    next if $i_evalue > $cutoff;
	    $domains{$domain}->{$sp}++;	    
	    push @{$domain_seqs{$domain}}, $gene if $dbh;
	}
    }
}
my @sp = sort keys %spnames;
print join("\t",
	   qw(DOMAIN TOTAL),
	   @sp),"\n";

for my $dom ( sort { $b->[1] <=> $a->[1] } 
	      map { [$_,sum(values %{$domains{$_}}), ] } 
	      keys %domains ) {
    if( $dbh ) {
	my $domname = $dom->[0];
	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$pfam_seqs_out/$domname.fa");
	my %seen;
	for my $id ( @{$domain_seqs{$domname}} ) {
	    next if $seen{$id}++;
	    if( my $seq = $dbh->get_Seq_by_acc($id) ) {
		$out->write_seq($seq);
	    } else {
		warn("cannot find $id in the DB\n");
	    }
	}
    }
    print join("\t", $dom->[0], $dom->[1],
	       map { $domains{$dom->[0]}->{$_} || 0 } @sp), "\n";    
}
