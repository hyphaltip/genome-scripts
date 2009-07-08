#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::Range;
my $version_out = 2;
my $debug = 0;
my $score = 10;
my $minlen = 50;
# assumes GFF3 input
GetOptions(
	   'gffv|gff:s' => \$version_out,
	   'v|debug!'   => \$debug,
	   's|score:s'  => \$score,
	   'm|minlen:i' => \$minlen,
	   );
if( $version_out == 3 ) {
	print "##gff-version 3\n";
}
my %genes;
my $gene;
while(<>) {
    next if /^\#/;
    chomp;
    my @row = split(/\t/,$_);
    my %group = map { split( /=/,$_) } split(/;/,pop @row);
    my $name;
    if( $row[2] eq 'gene' ) {
	push @{$genes{$gene->{'Name'}}}, $gene if $gene;
	$gene = undef;
    } elsif( $row[2] eq 'mRNA' ) {
	$name = $group{'ID'};	
	$gene = { 'Name' => $name,
		  'mRNA' => [@row],
		  'CDS'=> [],
	      };
    } elsif( $row[2] eq 'CDS' ) {
	push @{$gene->{'CDS'}}, [@row,$group{'Target'}];
    }
}
# fencepost
push @{$genes{$gene->{'Name'}}}, $gene if $gene;
my %unique;
for my $gene ( sort keys %genes ) {
    warn("gene is $gene\n") if $debug;
    my %chrom;
    for my $g ( @{$genes{$gene}} ) {
	my $seqid = $g->{'mRNA'}->[0];
	my $start = $g->{'mRNA'}->[3];
	my $end = $g->{'mRNA'}->[4];
	my $range = Bio::Range->new(-start => $start,
				    -end   => $end);
        # Technically this range includes introns, but let's not worry 
	# about it here, if there is a splice-site supported that is probably
	# okay information to include
	next if $range->length < $minlen;  
	push @{$chrom{$seqid}}, [$range,$g];
    }
    for my $chrom ( keys %chrom ) {
	my @sets = @{$chrom{$chrom}};
	
	if( @sets > 1 ) {
	    if( $debug ) {
		warn "BEFORE:\n";
		for my $s ( @sets ) {
		    my $s = $s->[1];
		    my $name = $s->{'Name'};
		    if( $version_out == 2 ) {
			warn join("\t", @{$s->{'mRNA'}}, 
				   sprintf("GenePrediction %s",
					   $name)),"\n";
			
		    } else {
			my @r = @{$s->{'mRNA'}};
			$r[2] = 'gene';
			warn join("\t",@r , 
				   sprintf("ID=%s.gene;Name=%s.gene",
					   $name,$name)),"\n";
			warn join("\t", @{$s->{'mRNA'}}, 
				   sprintf("ID=%s;Parent=%s.gene;Name=%s",
					   $name,$name,$name)),"\n";
		    }
		    for my $cds ( @{$s->{'CDS'}} ) {
			my @cds_i = @$cds;
			my $target = pop @cds_i;
			if( $version_out == 2) {
			    warn join("\t", @cds_i, 
				       sprintf("GenePrediction %s ; Target \"%s\"",
					       $name,$target)),"\n";
			} else {
			    warn join("\t", @cds_i, 
				       sprintf("Parent=%s;Target=%s",
					       $name,$target)),"\n";
			}
		    }
		}
	    }
	    # check for overlaps
	    # take longest one
	    my $changes = 1;
	    while( $changes ) {
		$changes = 0;
	      LOOP1: for( my $i = 0; $i < scalar @sets; $i++ ) {
		  my $s1 = $sets[$i]->[0];
		  for( my $j = $i+1; $j < scalar @sets; $j++ ) {
		      my $s2 = $sets[$j]->[0];
		      if( $s1->overlaps($s2) ) {
			  if( $s1->length > $s2->length ) {
			      splice(@sets,$j,1);
			  } else {
			      splice(@sets,$i,1);
			  }
			  $changes = 1;
			  last LOOP1;
		      }
		  }
	      }
	    }
	    for my $s ( @sets ) {
		$s->[1]->{'Name'} .= "-$chrom";
		if( $unique{$s->[1]->{'Name'}}++ ) {
		    $s->[1]->{'Name'} .= 
			sprintf(".copy%d",$unique{$s->[1]->{'Name'}});
		}
	    }
	} else {
	    $sets[0]->[1]->{'Name'} .= "-$chrom";
	}
	warn "AFTER:\n" if $debug;
	for my $s ( @sets ) {
	    $s = $s->[1];
	    my $name = $s->{'Name'};
	    if( $version_out == 2 ) {
		print join("\t", @{$s->{'mRNA'}}, sprintf("GenePrediction %s",
							  $name)),"\n";
	    } else {
		my @r = @{$s->{'mRNA'}};
		$r[2] = 'gene';
		print join("\t",@r , 
			   sprintf("ID=%s.gene;Name=%s.gene",
				   $name,$name)),"\n";
		print join("\t", @{$s->{'mRNA'}}, 
			   sprintf("ID=%s;Parent=%s.gene;Name=%s",
				   $name,$name,$name)),"\n";
	    }
	    for my $cds ( @{$s->{'CDS'}} ) {
		my $target = pop @$cds;
		if( $version_out == 2) {
		    print join("\t", @$cds, 
			       sprintf("GenePrediction %s ; Target \"%s\"",
				       $name,$target)),"\n";
		} else {
		    print join("\t", @$cds, 
			       sprintf("Parent=%s;Target=%s",
				       $name,$target)),"\n";
		}
	    }
	}
    }
}

