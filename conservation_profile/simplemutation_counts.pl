#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use File::Spec;
use Getopt::Long;
use List::Util qw(sum);
use Bio::SeqIO;
use constant GAP => '-';


# Algorithm
# -- process each alignment, 
#    count the number of different types of mutations

# mercator is 0-based
my $aln_file_name = 'output.mfa';
my $alnformat     = 'fasta';
my $alndir = 'alignments/pecan_alignments';
my $debug = 0;
my $ref_genome_index = 0;
my $outfile;
my $gaps_allowed = 0;
my $ambig_allowed = 0;

GetOptions(
	   'g|gaps!'          => \$gaps_allowed,
	   'f|alnfile:s'      => \$aln_file_name,
	   'af|format:s'      => \$alnformat,
	   'v|verbose!'       => \$debug,
	   'o|out|output:s'   => \$outfile,
	   'r|refidx|ref:i'   => \$ref_genome_index,
	   'd|a|i|input:s'    => \$alndir,
	   );

unless( -d $alndir || -r $alndir ) {
    die("cannot open alignment dir '$alndir'\n");
}


my %mutations;
my $total = 0;
opendir(DIR, $alndir) || die $!;
for my $aln_id ( sort { $a <=> $b } 
		 grep { -d "$alndir/$_" && /^\d+$/ } 
		 readdir(DIR) ) {
    warn("dir is $aln_id\n");
    my $fname = File::Spec->catfile($alndir,$aln_id,$aln_file_name);
    my $alnio = Bio::AlignIO->new(-file => $fname,
				  -format => $alnformat);
    if( my $aln = $alnio->next_aln ) {
	my $col_count    = $aln->length;

	# simplified data structure, we're going to use substr
	# for speed I hope
	my (@seqs,@ids);
	for my $seq ( $aln->each_seq ) {
	    push @seqs, $seq->seq;
	    push @ids, $seq->display_id; # just to double check
	    $seq = undef;
	}
	undef $aln;
	# process each column
	for( my $col_index = 0; $col_index < $col_count; $col_index++ ) {
	    my $ref_base = uc(substr($seqs[$ref_genome_index],$col_index,1));
	    if( $ref_base ne GAP && $ref_base ne 'N' ) { # we skip all ref-base gap columns
		# as these can't be written in ref coords
		# iterate thru each sequence
		my @alleles = ($ref_base);
		my $i = 0;		
		my ($score,$gap) = (0,0);
		my %alleles;
		$alleles{$ref_base}->[$i] = 1;
		my ($ambig,$gapped);
		for my $rest ( @seqs ) {
		    if( $i != $ref_genome_index ) { # that is NOT the ref
			my $base = uc(substr($seqs[$i],$col_index,1));
			$gapped++ if( $base eq GAP );
			$ambig++ if( $base eq 'N' );
			push @alleles, $base;
			$alleles{$base}->[$i] = 1;
		    }
		    $i++;			
		}
		# skip when there are gaps, unless they are allowed
		next if( $gapped && ! $gaps_allowed);
		next if( $ambig && ! $ambig_allowed);
		next if( @alleles != 3);
		$total++;
		
		if( scalar(keys %alleles) == 1 ) {
		    # identical
		    $mutations{'IDENTICAL'}++;
		} elsif( scalar (keys %alleles) == 3 ) {
		    # all different
		    $mutations{'ALL DIFFERENT'}++;
		} elsif( $alleles[0] eq $alleles[1] ) { # N.crassa=N.tetrasperma
		    # N.discreta mutation or 
		    # on N.crassa-N.tetrasperma ancestor
		    $mutations{'sp'}->{'N.DISCRETA'}->{$alleles[0]}->{$alleles[2]}++;
		} elsif( $alleles[0] eq $alleles[2] ) {
		    # N.tetrasperma mutation
		    $mutations{'sp'}->{'N.TETRASPERMA'}->{$alleles[0]}->{$alleles[1]}++;
		} elsif( $alleles[1] eq $alleles[2] ) {
		    # N.crass mutation
		    $mutations{'sp'}->{'N.CRASSA'}->{$alleles[1]}->{$alleles[0]}++;
		} else {
		    warn("unclassified @alleles\n");
		}		
	    }
	}
    } else {
	warn("No alignments in $fname\n");
    }
    warn("  done with ($aln_id)\n");   
    last if $debug;
}


printf "%s %d / %d (%5.2f%%)\n", 
    'IDENTICAL', $mutations{'IDENTICAL'}, $total, 
    100 * ($mutations{'IDENTICAL'} / $total);

printf "%s %d / %d (%5.2f%%)\n", 
    'ALL DIFFERENT', $mutations{'ALL DIFFERENT'}, $total, 
    100 * ($mutations{'ALL DIFFERENT'} / $total);

for my $sp ( sort keys %{$mutations{'sp'}} ) {
    my $total_cat = 0;
    my $d = $mutations{'sp'}->{$sp};
    
    my %all_bases;
    while( my ($ab,$db) = each %{$d} ) {	
	$all_bases{$ab}++;
	while( my ($ob, $ct) = each %{$db} ) {
	    $all_bases{$ob}++;
	    $total_cat += $ct;
	}
    }
    printf "%s %d / %d (%5.2f%%): \n", $sp,$total_cat, $total,
    100 * ($total_cat / $total);
    my @key_bases = sort keys %all_bases;
    for my $ab ( @key_bases ) {
	for my $ob ( @key_bases ) {
	    my $t = $d->{$ab}->{$ob} || 0;
	    printf "\t(anc) %s -> (obs) %s : %d (%5.2f%%)\n",
	    $ab, $ob, $t, 100 * ($t / $total_cat);
	}    
    }
    
    print "#\t",'ROW=ANCESTRAL COL=EXTANT',"\n";
    print join("\t", '',@key_bases, 'TOTAL'), "\n"; 
    my (%sums,$all);
    for(my $i = 0; $i < scalar @key_bases; $i++ ) {
	my @row;
	for(my $j = 0; $j < scalar @key_bases; $j++ ) {
	    push @row, $d->{$key_bases[$i]}->{$key_bases[$j]} || 0;
	    $sums{$key_bases[$j]} += $d->{$key_bases[$i]}->{$key_bases[$j]} || 0;
	}
	print join("\t", $key_bases[$i], @row, sum(@row)),"\n";
	$all += sum(@row);
    }
    print join("\t", '', (map { $sums{$_} || 0 } @key_bases), 
	       $all), "\n";
}

    
