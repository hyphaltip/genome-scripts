#!/usr/bin/perl -w
use strict;

=head1 NAME

gtf2gff3_3level_noCDS - convert gtf to GFF3 which is just exons and not CDS info

=head1 DESCRIPTION

This script is useful for Exons inferred from ESTs but no translation frame determined.

=head1 AUTHOR

Jason Stajich

=cut


my $debug = 0;
use Getopt::Long;
my ($transprefix,$prefix) = ( '','');
my $skipexons;
GetOptions('p|prefix:s' => \$prefix,
	   'tp|transprefix:s' => \$transprefix,
	   'v|debug!' => \$debug);
my %genes;
my %genes2alias;
my %transcript2name;
my %gene2name;
my $last_tid;
while(<>) {
    chomp;
    my $line = $_;
    my @line = split(/\t/,$_);
    next unless ($line[2] eq 'exon');
    $line[-1] =~ s/^\s+//;
    $line[-1] =~ s/\s+$//;
    my %set = map { split(/\s+/,$_,2) } split(/\s*;\s*/,pop @line);

    my ($gid,$tid,$tname,$gname,$alias) = 
	( map { $set{$_} =~ s/\"//g;
		$set{$_} }  
	  ( grep { exists $set{$_} } 
	    qw(gene_id transcript_id transcript_name gene_name aliases))
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
	if( ! @exons ) {
	    warn("no Exons\n");
	    next;
	}
	
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
	for my $exon ( @exons ) {
	    $exon = [$exon->[3],
		     join("\t", @$exon, sprintf("ID=%sexon%06d;Parent=%s",
						$prefix,
						$counts{'exon'}++,
						$mrna_id))];
	}
	print join("\n", ( map { $_->[1] } sort { $a->[0] <=> $b->[0] }
			   @exons)), "\n";
	
    }
}
