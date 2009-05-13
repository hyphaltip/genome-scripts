#!/usr/bin/perl -w

#  Turn BLAST M9 output into 
# GFF, presumably the file is
# EST-to-Genome mapping done by BLAT
# with the -out=blast9  option

use strict;
use Bio::DB::Fasta;
use Getopt::Long;
my $min_coverage = 0;
my $max_intron = 1000;
my $debug = 0;
GetOptions(
	   'mincoverage:i' => \$min_coverage,
	   'maxintron:i'   => \$max_intron,
	   'v|verbose!'    => \$debug,
	   );

my $dbdir = shift || die "need DBdir as first item!";
if ( ! -d $dbdir ) {
    die("need DBdir for first item not $dbdir\n");
}
my $estdb = Bio::DB::Fasta->new($dbdir);
my %groups;
my $group_SO_code = 'match';
my $hsp_SO_code   = 'HSP';
my $src = 'EST';

my $lastq;
while(<>) {
    next if /^\#/;
    my @line = split;

    my ($curq,$curh) = @line;
    if( $line[8] > $line[9] ) {
	push @line, -1;
	($line[8], $line[9]) = ( $line[9], $line[8]);
    } else { 
	push @line, 1;
    }
    if( defined $lastq && $lastq ne $curq ) {
	while(  my ($hit,$v) = each %groups ) {
	    my ($last_he, $last_strand, $last_qs, $last_hs,
		@transcripts,@hsps);
	    for my $exon ( sort { $a->[8] * $a->[12] <=> $b->[8] * $b->[12] } 
			   @$v ) {
		my ($q,$h, $pid, $hsplen, $gaps, $mismatches, 
		    $qs, $qe, $hs, $he, 
		    $evalue, $score, $strand) = @$exon;
		warn("#",join("\t", @$exon), "\n") if $debug;
		if( $last_he ) {
		    if( $last_strand != $strand ||
			$qs < $last_qs ) {
			# start new transcript
			push @transcripts, 
			[  $h, $src, $group_SO_code, 
			   $hs, $he, '.', $strand > 0 ? '+' : '-',
			   '.', $q,$qs,$qe];
		    } elsif( $strand > 0 ) {
			my $dist = ($hs - $last_he);
			if( $dist < 0 && $dist >= -3 ) {
			    warn("  hs was $hs\n") if $debug;
			    $hs += abs($dist)+1;
			    $dist = ($hs - $last_he);
			    warn("  hs is now $hs\n") if $debug;
			}
			if( $dist < -3 ) {
			    warn("dist is $dist, exons not in order for $q - $h [$strand] ($last_hs, $last_he, $hs) ($last_qs) -- $last_strand\n") if $debug;
			    next;
			} elsif( $dist > $max_intron ) {
			    push @transcripts, 
			    [  $h, $src, $group_SO_code, 
			       $hs, $he, '.', '+',
			       '.', $q,$qs,$qe];			    
			} elsif( $transcripts[-1]->[3] > $hs ||
				 $transcripts[-1]->[4] < $he ) {
			    if( $transcripts[-1]->[3] > $hs ) {
				$transcripts[-1]->[3] = $hs;
			    }
			    if( $transcripts[-1]->[4] < $he ) {
				$transcripts[-1]->[4] = $he;
			    }			  
			} else {
			    warn("weird, shouldn't be here\n");
			}
		    } else {
			my $dist = ($last_hs-$he);
			if( $dist < 0 && $dist > -3 ) {
			    warn("  he was $he\n") if $debug;
			    $he -= abs($dist)+1;
			    $dist = ($last_hs-$he);
			    warn("  he is now $he\n") if $debug;
			}
			if( $dist < 0 ) {
			    warn("exons not in order for $q - $h [$strand]($last_hs, $last_he, $last_qs\n) -- $last_strand") if $debug;
			    next;
			} elsif( $dist > $max_intron ) {
			    push @transcripts, 
			    [  $h, $src, $group_SO_code, 
			       $hs, $he, '.', '-','.', $q,$qs,$qe];
			} elsif( $transcripts[-1]->[3] > $hs ||
				 $transcripts[-1]->[4] < $he ) {
			    if( $transcripts[-1]->[3] > $hs ) {
				$transcripts[-1]->[3] = $hs;
			    }
			    if( $transcripts[-1]->[4] < $he ) {
				$transcripts[-1]->[4] = $he;
			    }
			    
			    if( $transcripts[-1]->[9] > $qs ){ 
				$transcripts[-1]->[9] = $qs;
			    }
			    if( $transcripts[-1]->[10] < $qe ) {
				$transcripts[-1]->[10] = $qe;
			    }
			} else {
			    warn("weird, shouldn't be here q=$q h=$h dist=$dist hs,he qs,qe ($hs,$he) ($qs,$qe)\n");
			}		    
		    }
		} else {
		    push @transcripts, [$h, $src, $group_SO_code, 
					$hs, $he, '.', 
					$strand > 0 ? '+' : '-', 
					'.', $q, $qs,$qe];
		}
		# preform the GFF string and include qs for
		# sorting when printing the HSPs in query order
		push @{$hsps[$#transcripts]},
		[ join("\t", $h, $src, $hsp_SO_code, $hs, $he, $score, 
		       $strand > 0 ? '+' : '-',
		       '.', sprintf("Target \"%s.%d %d %d\"", $q,
				    $#transcripts+1,$qs,$qe)), $qs];
		($last_he, $last_qs, $last_hs,$last_strand) = ( $he, $qs,
								$hs,$strand);
	    }
	    my $i= 0;
	    for my $transcript (@transcripts ) {
		
		my ($q,$qs,$qe) = splice(@$transcript,-3);
		my $len = $qe - $qs;
		my $estlen = $estdb->length($q);
		my $coverage = int(100 * ($len / $estlen));
		if( $coverage >= $min_coverage ) {
		    print join("\t",@$transcript,
			       sprintf("Target \"%s.%d %d %d\"; Length %d; Coverage %d; ESTLen %d",
				       $q,$i+1,$qs,$qe,
				       $len,
				       $coverage,
				       $estlen)),"\n";
		    # sort of query start (slot 1), 
		    # we already formed the GFF string in slot 0.
		    for my $hsp ( sort { $a->[1] <=> $b->[1] }
				  @{$hsps[$i]} ) {
			print $hsp->[0], "\n";
		    }
		}
		$i++;
	    }
	    warn("\n") if $debug;	    
	}	
	%groups = ();
#	last;
    }
    push @{$groups{$curh}}, [@line];
    $lastq = $curq;
}

