#!/usr/bin/perl -w
use strict;
my $fix = 0;
my $debug = 0;
use Getopt::Long;
my ($transprefix,$prefix) = ( '','');
my $skipexons;
GetOptions('fix!' => \$fix,
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

    my ($gid,$tid,$pid,$tname,$gname,$alias) = 
	( map { $set{$_} =~ s/\"//g;
		$set{$_} }  
	  ( grep { exists $set{$_} } 
	    qw(gene_id transcript_id protein_id
	       transcript_name gene_name aliases))
	  );    
    if( ! $tid ) {
	$tid = $last_tid;
    }
    if( defined $tid && $tid =~ /^\d+$/ ) { # JGI transcript ids are numbers only
	$tid = "t_$tid";
    }
    if( $tname ) {
	$transcript2name{$tid} = $tname;
    }
    if( ! $gid || ! $tid) {
	warn(join(" ", keys %set), "\n");
	die "Not GID or TID invalid GTF: $line \n";
    }
    if( $pid ) {
	$transcript2protein{$tid} = $pid;
    }
    # warn("tid=$tid pid=$pid gid=$gid tname=$tname gname=$gname\n");
    if( $fix ) {
	if( $tid =~ /(\S+)\.\d+$/) {
	    $gid = $1;
	}
    }
    if( $gname) {
	$gene2name{$gid} = $gname;
    }

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

print "##gff-version 3\n","##date-created ".localtime(),"\n";
my %counts;
for my $gid ( sort { $genes{$a}->{chrom} cmp $genes{$b}->{chrom} ||
			  $genes{$a}->{min} <=> $genes{$b}->{min}
		  } keys %genes ) {
    my $gene = $genes{$gid};
    my $gene_id = sprintf("%sgene%06d",$prefix,$counts{'gene'}++);
    my $aliases = $genes2alias{$gid};
    my $gname   = $gene2name{$gid};
    if( $gname ) {
	if( $aliases) {
	    $aliases = join(",",$gid,$aliases);	
	} else {
	    $aliases = $gid;
	}
    } else {
	$gname = $gid;
    }
    $gname = sprintf('"%s"',$gname) if $gname =~ /[;\s,]/;
    my $ninth  = sprintf("ID=%s;Name=%s",$gene_id, $gname);
    if( $aliases ) {
	$ninth .= sprintf(";Alias=%s",$aliases);
    }
    print join("\t", ( $gene->{chrom}, 
		       $gene->{src},
		       'gene',
		       $gene->{min},
		       $gene->{max},
		       '.',
		       $gene->{strand},
		       '.',
		       $ninth)),"\n";
    while( my ($transcript,$exons) = each %{$gene->{'transcripts'}} ) {
	my $mrna_id = sprintf("%smRNA%06d",$prefix,$counts{'mRNA'}++);
	my @exons = grep { $_->[2] eq 'exon' } @$exons;
	my @cds   = grep { $_->[2] eq 'CDS'  } @$exons;	
	if( ! @cds ) {
	    warn("no CDS found in $mrna_id ($gid,$transcript) ", join(",", map { $_->[2] } @$exons),"\n") if $debug;
	    next;
	}
	my ($start_codon)   = grep { $_->[2] eq 'start_codon'  } @$exons;
	my ($stop_codon)   = grep { $_->[2] eq 'stop_codon'  } @$exons;
	
	if( ! @exons ) {
	  for my $e ( @cds ) {
	    push @exons, [@$e];
	    $exons[-1]->[2] = 'exon';
	  }
	}
	
	my $proteinid = $transcript2protein{$transcript};
	my ($chrom,$src,$strand,$min,$max);
	for my $exon ( @exons ) {
	    $chrom = $exon->[0] unless defined $chrom;
	    $src   = $exon->[1] unless defined $src;
	    $min   = $exon->[3] if( ! defined $min || $min > $exon->[3]);
	    $max   = $exon->[4] if( ! defined $max || $max < $exon->[4]);
	    $strand = $exon->[6] unless defined $strand;	
	}

	my $strand_val = $strand eq '-' ? -1 : 1;
	my $transname = $transprefix.$transcript;
	my $transaliases = $transcript;
	if( exists $transcript2name{$transcript} ) {
	    $transname = $transcript2name{$transcript};
	    $transname = sprintf('"%s"',$transname) if $transname =~ /[;\s,]/;
	}
	my $mrna_ninth = sprintf("ID=%s;Parent=%s;Name=%s",
				 $mrna_id,$gene_id,$transname);
	if( $transaliases && $transaliases ne $transname ) {
	    $mrna_ninth .= sprintf(";Alias=%s",$transaliases);
	}
	print join("\t",($chrom,
			 $src,
			 'mRNA',
			 $min,
			 $max,
			 '.',
			 $strand,
			 '.',
			 $mrna_ninth,
			 )),"\n";
	
	my ($translation_start, $translation_stop);
	
	# order 5' -> 3' by multiplying start by strand

	@cds = sort { $a->[3] * $strand_val <=>
		      $b->[3] * $strand_val } @cds;

	if ( $debug) {
	  warn("CDS order is :\n");
	  for my $c ( @cds ) {
	    warn(join("\t", @$c), "\n");
	  }
	}
	
	if( $stop_codon ) {
	    if( $strand_val > 0 ) {		
		warn("stop codon is ", join("\t", @{$stop_codon}), "\n") if $debug;
		$cds[-1]->[4] = $stop_codon->[4];
		$translation_stop = $stop_codon->[4];
	    } else {
		warn("stop codon is ", join("\t", @{$stop_codon}), "\n") if $debug;
		$cds[-1]->[3] = $stop_codon->[3];
		$translation_stop = $stop_codon->[3];
	    }
	  } else {
	    $translation_stop = ($strand_val > 0) ? $cds[-1]->[4] : $cds[-1]->[3];
	  }
	if( $start_codon ) {
	    if( $strand_val > 0 ) {		
		warn("start codon is ", join("\t", @{$start_codon}), "\n") if $debug;
		$cds[0]->[3] = $start_codon->[3];
		$translation_start = $start_codon->[3];
	    } else {
		warn("start codon is ", join("\t", @{$start_codon}), "\n") if $debug;
		$cds[0]->[4] = $start_codon->[4];
		$translation_start = $start_codon->[4];
	    }
	  } else {
	    $translation_start = ($strand_val > 0) ? $cds[0]->[3] : $cds[0]->[4];
	  }
	if ( $debug) {
	  warn("CDS order is after :\n");
	  for my $c ( @cds ) {
	    warn(join("\t", @$c), "\n");
	  }
	}

	for my $cds_i ( @cds ) {
	  my $exon_ninth = sprintf("ID=%scds%06d;Parent=%s",
				   $prefix,
				   $counts{'CDS'}++,
				   $mrna_id);
	  if ( $proteinid ) {
	    $proteinid = sprintf('"%s"',$proteinid) if $proteinid =~ /[;\s,]/;
	    $exon_ninth .= sprintf(";Name=%s",$proteinid);
	  }
	  $cds_i = [$cds_i->[3], join("\t", @$cds_i, $exon_ninth)];
	}
	if ( $debug) {
	  warn("CDS order after ninth column is :\n");
	  for my $c ( @cds ) {
	    warn(join("\t", @$c), "\n");
	  }
	}

	# making some utrs
	my %utrs;
	for my $exon ( sort { $a->[3] * $strand_val <=> 
			      $b->[3] * $strand_val }
		       @exons ) {
	  # how many levels deep can you think? ...
	  if ( defined $translation_start && defined $translation_stop ) {
	    if ( $strand_val > 0 ) {
	      # 5' UTR on +1 strand
	      if ( $translation_start > $exon->[3] ) {
		if ( $translation_start > $exon->[4] ) {
		  # whole exon is a UTR so push it all on
		  push @{$utrs{'5utr'}},
		    [ $exon->[3],
		      join("\t",
			   ( $exon->[0],
			     $exon->[1],
			     'five_prime_utr',
			     $exon->[3],
			     $exon->[4],
			     '.',
			     $strand,
			     '.',
			     sprintf("ID=%sutr5%06d;Parent=%s",
				     $prefix,
				     $counts{'5utr'}++,
				     $mrna_id)))];
		} else {
		  # push the partial exon up to the start codon
		  push @{$utrs{'5utr'}}, 
		    [ $exon->[3],
		      join("\t",
			   $exon->[0],
			   $exon->[1],
			   'five_prime_utr',
			   $exon->[3],
			   $translation_start - 1,
			   '.',
			   $strand,
			   '.',
			   sprintf("ID=%sutr5%06d;Parent=%s",
				   $prefix,
				   $counts{'5utr'}++,
				   $mrna_id)
			  )];
		}
	      }
	      #3' UTR on +1 strand
	      if ( $translation_stop < $exon->[4] ) {
		if ( $translation_stop < $exon->[3] ) {
		  # whole exon is 3' UTR
		  push @{$utrs{'3utr'}},
		    [ $exon->[3],
		      join("\t",
			   ( $exon->[0],
			     $exon->[1],
			     'three_prime_utr',
			     $exon->[3],
			     $exon->[4],
			     '.',
			     $strand,
			     '.',
			     sprintf("ID=%sutr3%06d;Parent=%s",
				     $prefix,
				     $counts{'3utr'}++,
				     $mrna_id)))];
		} else { 
		  # make UTR from partial exon
		  push @{$utrs{'3utr'}},
		    [ $exon->[3],
		      join("\t",
			   ( $exon->[0],
			     $exon->[1],
			     'three_prime_utr',
			     $stop_codon +1,
			     $exon->[4],
			     '.',
			     $strand,
			     '.',
			     sprintf("ID=%sutr3%06d;Parent=%s",
				     $prefix,
				     $counts{'3utr'}++,
				     $mrna_id)))];
		}
	      } 
	    } else {
	      # 5' UTR on -1 strand
	      if ( $translation_start < $exon->[4] ) {
		if ( $translation_start < $exon->[3] ) {
		  # whole exon is UTR
		  push @{$utrs{'5utr'}},
		    [ $exon->[3],
		      join("\t",
			   $exon->[0],
			   $exon->[1],
			   'five_prime_utr',
			   $exon->[3],
			   $exon->[4],
			   '.',
			   $strand,
			   '.',
			   sprintf("ID=%sutr5%06d;Parent=%s",
				   $prefix,
				   $counts{'5utr'}++,
				   $mrna_id)) ];
		} else {
		  # push on part of exon up to the start codon 
		  push @{$utrs{'5utr'}}, 
		    [ $exon->[3],
		      join("\t",$exon->[0],
			   $exon->[1],
			   'five_prime_utr',
			   $translation_start + 1,
			   $exon->[4],
			   '.',
			   $strand,
			   '.',
			   sprintf("ID=%sutr5%06d;Parent=%s",
				   $prefix,
				   $counts{'5utr'}++,
				   $mrna_id))];
		}
	      }		
	      #3' UTR on -1 strand
	      if ( $translation_stop > $exon->[3] ) {
		if ( $translation_stop > $exon->[4] ) {
		  # whole exon is 3' UTR
		  push @{$utrs{'3utr'}},
		    [ $exon->[3],
		      join("\t",
			   ( $exon->[0],
			     $exon->[1],
			     'three_prime_utr',
			     $exon->[3],
			     $exon->[4],
			     '.',
			     $strand,
			     '.',
			     sprintf("ID=%sutr3%06d;Parent=%s",
				     $prefix,
				     $counts{'3utr'}++,
				     $mrna_id)))];
		} else { 
		  # make UTR from partial exon
		  push @{$utrs{'3utr'}},
		    [ $exon->[3],
		      join("\t",
			   ( $exon->[0],
			     $exon->[1],
			     'three_prime_utr',
			     $exon->[3],
			     $translation_stop -1,
			     '.',
			     $strand,
			     '.',
			     sprintf("ID=%sutr3%06d;Parent=%s",
				     $prefix,
				     $counts{'3utr'}++,
				     $mrna_id)))];
		}
	      }
	    }
	  }
	  $exon = [$exon->[3],
		   join("\t", @$exon, sprintf("ID=%sexon%06d;Parent=%s",
					      $prefix,
					      $counts{'exon'}++,
					      $mrna_id))];
	}
	if( $strand_val > 0 ) {
	    if( exists $utrs{'5utr'} ) {
		
		print join("\n", map { $_->[1] } sort { $a->[0] <=> $b->[0] }
			   @{$utrs{'5utr'}}), "\n";
	    }
	} else {
	    if( exists $utrs{'3utr'} ) {
		print join("\n", map { $_->[1] } sort { $a->[0] <=> $b->[0] }
			   @{$utrs{'3utr'}}), "\n";
	    }
	}
	
	print join("\n", ( map { $_->[1] } sort { $a->[0] <=> $b->[0] }
			   @exons, @cds)), "\n";
	if( $strand_val > 0 ) {
	  if( exists $utrs{'3utr'} ) {
	    print join("\n", map { $_->[1] } sort { $a->[0] <=> $b->[0] }
		       @{$utrs{'3utr'}}), "\n";
	  }
	} else {
	  if( exists $utrs{'5utr'} ) {
	    print join("\n", map { $_->[1] } sort { $a->[0] <=> $b->[0] }
		       @{$utrs{'5utr'}}), "\n";
	    }
	}
    }
}
