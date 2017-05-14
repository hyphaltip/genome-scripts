#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqIO;

=head1 NAME 

pfamcopy_compare - compare Pfam copy number counts for a collection of species

=head1 USAGE

pfamcopy_compare.pl pfam_outputdir

=head1 DESCRIPTION

1. hmmpfam (v3) to table of counts per species
2. If seqs provided will make a file for each domain with seq domain cut out

=head1 AUTHOR

Jason Stajich E<lt>jason.stajich[at]ucr.eduE<gt>

=cut

my $ext = '\.domtbl\.tab';
my $cutoff = 1e-4;
my $pfam_seqs_out = 'pfam_seq_out';
my $dbdir;
my $percent;
GetOptions('c|cutoff:s' => \$cutoff,
	   'e|ext:s'    => \$ext,
	   'db:s' => \$dbdir,
	   'p|percentage!' => \$percent,
	   'o|out:s'=> \$pfam_seqs_out);

my $dir = shift || 'pfam';

my $dbh;
if( $dbdir && (-f $dbdir || -d $dbdir) ) {
    $dbh = Bio::DB::Fasta->new($dbdir);
    mkdir($pfam_seqs_out);
}
opendir(DIR,$dir) || die "$dir: $!";
my (%domains,%spnames,%domain_seqs,%domains_genes,%spgenes,%genes);
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

	    if( $percent ) {
		$domains{$domain}->{$sp}->{$gene}++;
	    } else {
		$domains{$domain}->{$sp}++;	
	    }

	    $genes{$sp}->{$gene}++;
	    push @{$domain_seqs{$domain}}, [$gene,$env_start,$env_end] if $dbh;
	}
    }
}
my @sp = sort keys %spnames;
for my $sp ( @sp ) {
    $spgenes{$sp} = scalar keys %{$genes{$sp}};
}

if( $percent ) {
    print join("\t",
	       qw(DOMAIN),
	       @sp),"\n";
} else {
    print join("\t",
	       qw(DOMAIN TOTAL),
	       @sp),"\n";
}

my @dom_order;
if( $percent ) {
    @dom_order = map { [$_,$_] } sort keys %domains;
} else {
    @dom_order = sort { $b->[1] <=> $a->[1] } 
    map {  [$_,sum(values %{$domains{$_}}), ] }
    keys %domains;
}	     
for my $dom ( @dom_order ) {
    if( $dbh && ! $percent ) {
	my $domname = $dom->[0];
	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => ">$pfam_seqs_out/$domname.fa");
	my %seen;
	for my $dat ( sort { $a->[0] cmp $b->[0] ||
				 $a->[1] <=> $b->[1] ||
				 $a->[2] <=> $b->[2] } 
		      @{$domain_seqs{$domname}} ) {
	    my ($id,$start,$end) = @$dat;	    
	    my $domain_n = ++$seen{$id};
	    if( my $seq = $dbh->seq($id, $start => $end) ) { 
		my $domainseq = Bio::PrimarySeq->new
		    (-seq => $seq, 
		     -id  => sprintf("%s_D%02d",$id,$domain_n),
		     -desc => sprintf("%s %d..%d",$id, $start,$end));
		$out->write_seq($domainseq);
	    } else {
		warn("cannot find $id in the DB\n");
	    }
	}
    }
    if( $percent ) {
	print join("\t", $dom->[0], 
		   map {my $ct = scalar keys %{$domains{$dom->[0]}->{$_} || {}};
			sprintf("%.2f", 100 * ($ct / $spgenes{$_}));
		   } @sp), "\n";

    } else {
	print join("\t", $dom->[0], $dom->[1],
		   map { $domains{$dom->[0]}->{$_} || 0 } @sp), "\n";    
    }
}
