#!/usr/bin/perl -w
use strict;
# author jason stajich <stajich@berkeley.edu>
# description:
# this script will turn GFF (2 or 3) into gtf (gff 2.5) suitable for twinscan training 
# by adding start and stop codons (removing 1st and last codons from the annotated CDS).
# There are probably a few problems with this approach, but it works well enough and I am assuming the
# validate_gtf.pl script from the Brent lab catches any problems to weed out.
#
# It currently expects GFF2 data to have either Transcript or GenePrediction group tag
use Bio::DB::Fasta;
use Getopt::Long;

my $debug = 0;
GetOptions (
	    'v|version!' => \$debug);

my $db= Bio::DB::Fasta->new(shift @ARGV);
# Frame is calculated as (3 - ((length-frame) mod 3)) mod 3
my @order;
my %gene;
my %seen;
while(<>) {
    my @line = split(/\t/,$_);
    next if uc($line[2]) ne 'CDS';
    my $last = pop @line;
    chomp($last);
    my $group;
    if( $last =~ /(Transcript|GenePrediction)\s+(\S+)/ ) {
	($group) = $2;
    } elsif( $last =~ /Parent=([^;]+);?/) {
	$group = $1;
	$group =~ s/Model\.//;
    } 
     
    if( ! $seen{$group}++ ) {
	push @order, $group;
    }
    push @{$gene{$group}}, [ @line, 
			     sprintf('gene_id "%s"; transcript_id "%s";',
				     $group, "$group.1")];
}
for my $gene ( @order ) {
    my @ordered_cds = ( map { $_->[1] }
			sort { $a->[0] <=> $b->[0]}
			map { [$_->[3] * ($_->[6] eq '-' ? -1 : 1), $_] }
			@{$gene{$gene}} );
    my ($fexon,$lexon) = ($ordered_cds[0], $ordered_cds[-1]);
    if( $fexon->[6] eq '+' ) {
	print join("\t", $fexon->[0], $fexon->[1], 'start_codon',
		   $fexon->[3],
		   $fexon->[3] + 2,
		   $fexon->[5],
		   $fexon->[6],
		   $fexon->[7],
		   $fexon->[8]), "\n";
	$fexon->[3] += 3;
    } else {
	print join("\t", $fexon->[0], $fexon->[1], 'start_codon',
		   $fexon->[4]-2,
		   $fexon->[4],
		   $fexon->[5],
		   $fexon->[6],
		   $fexon->[7],
		   $fexon->[8]), "\n";	
	$fexon->[4] -= 3;
    }
    
    if($lexon->[6] eq '-' ){
	if( $debug ) {
	    my $last_codon = $db->seq($lexon->[0], 
				      $lexon->[3]+ 2 => $lexon->[3] );
	    
	    my $next_last_codon = $db->seq($lexon->[0], $lexon->[3]+5 => $lexon->[3]+3);
	    print "last_codon $last_codon : next_last_codon $next_last_codon\n";
	    print "last exon ",$db->seq($lexon->[0],
					$lexon->[4] => $lexon->[3]),"\n"; 
	}
	push @ordered_cds, [$lexon->[0], $lexon->[1], 'stop_codon',
			    $lexon->[3],
			    $lexon->[3] + 2,
			    $lexon->[5],
			    $lexon->[6],
			    $lexon->[7],
			    $lexon->[8],
			    ];
	$lexon->[3] += 3;
	if( $debug ) {
	    print "last exon ",$db->seq($lexon->[0],
					$lexon->[4] => $lexon->[3]),"\n"; 
	}
    } else {
	if( $debug ) {
	    my $last_codon = $db->seq($lexon->[0], $lexon->[4]-2 => $lexon->[4]);
	    my $next_last_codon = $db->seq($lexon->[0], $lexon->[4]-5 => $lexon->[4]-3);
	    print "last_codon $last_codon : next_last_codon $next_last_codon\n";
	    print "last exon ",$db->seq($lexon->[0],
					$lexon->[3] => $lexon->[4]),"\n"; 
	}
	push @ordered_cds,[$lexon->[0], $lexon->[1], 'stop_codon',
		   $lexon->[4]-2,
		   $lexon->[4],
		   $lexon->[5],
		   $lexon->[6],
		   $lexon->[7],
		   $lexon->[8]];
	$lexon->[4] -= 3;
	if( $debug )  {
	    print "last exon ",$db->seq($lexon->[0],
					$lexon->[3] => $lexon->[4]),"\n"; 
	}
    }
    for my $cds ( @ordered_cds ) {
	print join("\t", @$cds), "\n";
    }


}
