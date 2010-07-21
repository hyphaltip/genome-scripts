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
my %gene2name;
my %transcript2protein;
my $last_tid;
while(<>) {
    chomp;
    my $line = $_;
    my @line = split(/\t/,$_);
    next unless ($line[2] eq 'exon' );
    $line[-1] =~ s/^\s+//;
    $line[-1] =~ s/\s+$//;
    my %set = map { split(/\s+/,$_,2) } split(/\s*;\s*/,pop @line);

    my ($gid,$tid,$alias,$code,$exonnum,$tss_id) = 
	( map { defined $set{$_} && $set{$_} =~ s/\"//g;
		$set{$_} }  
	  qw(gene_id transcript_id nearest_ref class_code exon_number tss_id)
	  );    
    if( ! $tid ) {
	$tid = $last_tid;
    }
    if( defined $tid && $tid =~ /^\d+$/ ) { # JGI transcript ids are numbers only
	$tid = "t_$tid";
    }
    
    if( ! $gid || ! $tid) {
	warn(join(" ", keys %set), "\n");
	die "Not GID or TID invalid GTF: $line \n";
    }
    # warn("tid=$tid pid=$pid gid=$gid tname=$tname gname=$gname\n");
    if( $fix ) {
	if( $tid =~ /(\S+)\.\d+$/) {
	    $gid = $1;
	}
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
    push @{$genes{$gid}->{transcripts}->{$tid}}, [@line,
						  {'exon_number' => $exonnum,
						   'exon_class'  => $code,
						   }];
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
	
	# order 5' -> 3' by multiplying start by strand
	@exons = sort { $a->[3] * $strand_val <=>
			    $b->[3] * $strand_val } @exons;		
	for my $e ( @exons ) {
	    my $last = pop @$e;
	    print join("\t", @$e, sprintf("ID=exon%07d;Parent=%s;exon_num=%d;class=\"%s\"",
					  $counts{'exon'}++, 
					  $mrna_id, 
					  $last->{exon_number},
					  $last->{exon_class})),"\n";
	}
    }
    last if $debug;
}
