#!/usr/bin/perl -w
use strict;

my $fix = 0;
my $debug = 0;
use Getopt::Long;
my ($transprefix,$prefix) = ( '','');
my $skipexons;
GetOptions('fix!' => \$fix, # get point name
	   's|skipexons:s' => \$skipexons,
	   'p|prefix:s' => \$prefix,
	   'tp|transprefix:s' => \$transprefix,
	   'v|debug!' => \$debug);
my %genes;
my %genes2alias;
my %transcript2name;
my %gene2name;
my %transcript2protein;
my $last_tid;
while(<>) {
    chomp;
    my $line = $_;
    my @line = split(/\t/,$_);
    next unless ($line[2] eq 'CDS'  || 
		 $line[2] eq 'exon' || 
		 $line[2] eq 'stop_codon' ||
		 $line[2] eq 'start_codon'
		 );
    $line[-1] =~ s/^\s+//;
    $line[-1] =~ s/\s+$//;
    my %set = map { split(/\s+/,$_,2) } split(/\s*;\s*/,pop @line);

    my ($gid,$tid,$pid,$tname,$gname,$alias,$name,$transID,$protID) = 
	( map { 
	    my $rc = undef;
	    if( exists $set{$_} ) {
		$set{$_} =~ s/\"//g;
		$rc = $set{$_}; 
	    } 
	    $rc;
	} qw(gene_id transcript_id protein_id 
	     transcript_name gene_name aliases name transcriptId proteinId));
	      

    $tid = $transID if ! $tid && $transID;
    $pid = $protID if ! $pid && $protID;
    $gname = $name if ! $gname && $name;
    $gname = $tname if ! $gname && $tname;
        
    $gid ||= $gname;
    $tid ||= $gid.".T";

    if( ! defined $genes{$gid}->{min} ||
	$genes{$gid}->{min} > $line[3] ) {
	$genes{$gid}->{min} = $line[3];
    }
    if( ! defined $genes{$gid}->{max} ||
	$genes{$gid}->{max} < $line[4] ) {
	$genes{$gid}->{max} = $line[4];
    }
    if( ! defined $genes{$gid}->{strand} ) {
	$genes{$gid}->{strand} = $line[6];
    }
    if( ! defined $genes{$gid}->{chrom} ) {
	$genes{$gid}->{chrom} = $line[0];
    }
    if( ! defined $genes{$gid}->{src} ) {
	$genes{$gid}->{src} = $line[1];
    }
    if( defined $alias ) {
	$genes2alias{$gid} = join(',',split(/\s+/,$alias));
    }
    push @{$genes{$gid}->{transcripts}->{$tid}}, [@line];
    $last_tid = $tid;    
}

my %counts;
for my $gid ( sort { $genes{$a}->{chrom} cmp $genes{$b}->{chrom} ||
			  $genes{$a}->{min} <=> $genes{$b}->{min}
		  } keys %genes ) {
    my $gene = $genes{$gid};
    while( my ($transcript,$exonsref) = each %{$gene->{'transcripts'}} ) {
	my $cds_count = scalar grep { $_->[2] eq 'exon' } @$exonsref;
	my %ct = ('exon' => 1, 'CDS' => 1);
	for my $ex ( sort { ( $a->[3] * ($a->[6] eq '-' ? -1 : 1) <=> 
			      $b->[3] * ($b->[6] eq '-' ? -1 : 1)) 
			} @$exonsref ) {
	    my $ex_type;
	    if( $ex->[2] =~ /CDS|exon/) {
		if( $cds_count == 1 ) {
		    $ex_type = 'single';
		} elsif( $ct{$ex->[2]} == 1 ) {
		    $ex_type = 'initial';
		} elsif( $ct{$ex->[2]} == $cds_count ) {
		    $ex_type = 'terminal';
		} else { 
		    $ex_type = 'internal';
		}
		$ct{$ex->[2]}++;
	    }

	    if( $ex->[2] eq 'exon' ) {
		$ex->[8] = sprintf("exontype \"%s\"; transcript_id \"%s\"; gene_id \"%s\";",
				   $ex_type,$transcript, $gid);
	    } elsif( $ex->[2] eq 'CDS' ) {
		$ex->[8] = sprintf("exontype \"%s\"; transcript_id \"%s\"; gene_id \"%s\";",
				   $ex_type,$transcript, $gid);
	    } else { 
		# start or stop codon
		$ex->[8] = sprintf("transcript_id \"%s\"; gene_id \"%s\";",
				   $transcript, $gid);
	    }
       
	    print join("\t", @$ex), "\n";	    
	}
    }
}
