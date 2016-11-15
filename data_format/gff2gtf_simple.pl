#!/usr/bin/perl -w
use strict;
# author jason stajich <stajich@berkeley.edu>
# description:
# this script will turn GFF (2 or 3) into gtf (gff 2.5) suitable for htseq
use Getopt::Long;

my $debug = 0;
my $show_exontype = 0;
GetOptions (
    'v|verbose!' => \$debug,
    'show|exon!' => \$show_exontype,
    );

my $db;
if( $debug ) {
    $db = Bio::DB::Fasta->new(shift @ARGV);
}
# Frame is calculated as (3 - ((length-frame) mod 3)) mod 3
my @order;
my %gene;
my %seen;
my %id2name;
while(<>) {
    next if /^#/;
    my @line = split(/\t/,$_);
    my $last = pop @line;
    chomp($last);
    my $group;    
    if( uc($line[2]) eq 'MRNA' ||
	uc($line[2]) eq 'GENE') {
	my ($id,$name);
	if( $last =~ /ID=([^;]+);?/ ) {
	    $id = $1;
	}
	if( $last =~ /Name=([^;]+);?/ ) {
	    $name = $1;
	}
	#warn("id is $id\n");
	if( defined $id && defined $name ) {
	    $id2name{uc $line[2]}->{$id} = $name; # map the id (typically just a numeric) to the name (string)
	    if( uc($line[2]) eq 'MRNA' && 
		$last =~ /Parent=([^;]+)/ ) {
		# map the mRNA ID to the Gene (numeric) ID
		$id2name{MRNA_PARENT}->{$id} = $1;
	    }
	}
#	warn($_);
    }
    next unless uc $line[2] eq 'CDS' || $line[2] eq 'exon';
    if( $last =~ /(ID|Name|Transcript|GenePrediction)\s+(\S+)/ ) {
	($group) = $2;
    } elsif( $last =~ /Parent=([^;]+);?/) { # gff3
	$group = $1;
	$group =~ s/Model\.//;
    } 
    if( ! $group ) {
	warn("no group in $_\n");
	next;
    } 
    my $tid = "$group.1";
    my $gid = $group;    
    if( exists $id2name{'MRNA'}->{$group} ) {
	$tid = $id2name{'MRNA'}->{$group};
    } else {
	warn "cannot find MRNA group $group\n";
	die;
    }
	
    if( exists $id2name{'MRNA_PARENT'}->{$group} ) {
	my $geneid = $id2name{'MRNA_PARENT'}->{$group};
	if( exists $id2name{'GENE'}->{$geneid} ) {
	    $gid = $id2name{'GENE'}->{$geneid};
	}
    } else {
	warn("cannot find gene group for $group\n");
	die;
    }
    if( ! $seen{$gid}++ && uc($line[2]) eq 'CDS') {
	push @order, $gid;
	warn("adding $gid\n");
    }
	# warn("adding $gid\n");
    push @{$gene{$gid}}, [ @line, 
			   sprintf('gene_id "%s"; transcript_id "%s";',
#			       'transcript_id "%s"; gene_id "%s";',
				   $gid, $tid)];
}
if( ! @order ) {
	@order = sort keys %gene;
}
for my $gene ( @order ) {
	warn("gene is $gene\n");
    my @ordered_cds = ( map { $_->[1] }
			sort { $a->[0] <=> $b->[0]}
			map { [$_->[3] * ($_->[6] eq '-' ? -1 : 1), $_] }
			@{$gene{$gene}} );
    my $i = 1;
    my $count = scalar @ordered_cds;
    for my $cds ( @ordered_cds ) {
	print join("\t", @$cds), "\n";
    }
}
