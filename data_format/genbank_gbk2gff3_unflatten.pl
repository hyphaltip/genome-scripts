#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Bio::SeqFeature::Tools::Unflattener;

# generate an Unflattener object
my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;
$unflattener->error_threshold(1);

use constant MRNA => 'mRNA';
use constant GENE => 'gene';
use constant CDS  => 'cds';
use constant EXON  => 'exon';
use constant MIN_LENGTH => 10_000;

my $SRC = "genbank";

my $dir = "raw/genbank";
my $odir = 'genomes';
opendir(DIR, $dir) || die $!;


for my $subdir ( readdir(DIR) ) {    
    next unless -d "$dir/$subdir";
    opendir(SUBDIR,"$dir/$subdir");
    my $cdsnum = 1;
    my $first = 1;
    my (@ALL,$outfh,@seqs);
    for my $file ( readdir(SUBDIR)) {
	next unless $file =~ /(\S+)\.gb[sk]$/;
	warn("$file\n");
	my ($species,$pref) = &get_species_gbkfile("$dir/$subdir/$file");
	if( $first ) {
	    mkdir("$odir/$species") unless -d "$odir/$species";
	    next if( -f "$odir/$species/$species.gff3" && 
		     ! -z "$odir/$species/$species.gff3");
	    next if (-f "$odir/$species/$species.gff3.gz");
	    
	    warn("species is $species\n");
	    open($outfh, ">$odir/$species/$species.gff3") || die $!;
	    print $outfh "##gff-version 3\n";
	    print $outfh "#date ". localtime(time)."\n";
	}
	$first = 0;

	my $seqio = Bio::SeqIO->new(-format => 'genbank',
				    -file   => "$dir/$subdir/$file");

	my %genes;

	while( my $seq = $seqio->next_seq ) {
	    my $desc = $seq->description;
	    if( $desc =~ /(?:chromosome|contig)\s+([\w\d\.]+)/) {
		$seq->display_id("$pref\_$1");
	    }

	    my @top_sfs = $seq->get_SeqFeatures;
	    unless( $seq->length > MIN_LENGTH &&
		    grep { $_->primary_tag eq 'gene' } @top_sfs) {
		next;
	    }
	    warn("desc is $desc id is ", $seq->display_id,"\n");
	    # push @seqs, $seq;
	    print $outfh  join("\t",
			       $seq->display_id,
			       'chromosome',
			       'scaffold',
			       1, $seq->length,
			       '.','.','.',
			       sprintf("ID=%s;Name=%s;Accession=%s.%d",
				       $seq->display_id, $seq->display_id, 
				       $seq->accession_number, 
				       $seq->seq_version)),"\n";

	    # get top level unflattended SeqFeatureI objects
	    $unflattener->unflatten_seq(-seq=>$seq,
#					-group_tag=>'locus_tag',
					-use_magic=>1);
	    # $out->write_seq($seq);
	    my $i = 0;

	    foreach my $f (@top_sfs) {
		next unless $f->primary_tag eq 'gene';
		next if $f->has_tag('pseudo');
		my $primarytag = $f->primary_tag;
		next unless( $primarytag eq 'CDS' || 
			     $primarytag eq 'gene' || 
			     $primarytag eq 'mRNA');
		my $genestr;
		my ($min,$max,$strand) = ($f->start,$f->end, $f->strand);
		my ($pname,%genexrefs, @genexrefs_a);

		if( $f->has_tag('db_xref') ) { 
		    for my $xref ( $f->get_tag_values('db_xref') ) {
			my ($xref_src,$xref_id) = split(':',$xref);
			$genexrefs{$xref_src} = $xref_id;
			push @genexrefs_a, &escape($xref);
		    } 
		}
		for my $ptag ( qw(locus_tag gene name) ) {
		    if( $f->has_tag($ptag) ) {
			($pname) = $f->get_tag_values($ptag);
			last;
		    }
		}
		$pname = $genexrefs{GeneID} if( exists $genexrefs{'GeneID'} && 
						! defined $pname);

		unless( defined $pname ) {
		    warn("cannot find pname in ", $f->gff_string, "\n");
		    last;
		}

		my %refs;
		push @ALL, [$seq->display_id,
			    $SRC,
			    GENE, 
			    $f->start,
			    $f->end,
			    '.',
			    $f->strand < 0 ? '-' : '+',
			    '.',
			    { 'ID' => "$pname.gene", 'Name'=>$pname} ];

		if( @genexrefs_a ) {
		    $ALL[-1]->[8]->{'Dbxref'} = join(",", @genexrefs_a);
		}

		my $mrnact = 1;
		my @mrnas = $f->get_SeqFeatures;
		for my $mRNA ( @mrnas ) {
		    my $mrna_name = $pname;
		    if( @mrnas > 1 ) {
			$mrna_name = "$pname.$mrnact"; 
		    }
		    my $x = [$seq->display_id,
			     $SRC,
			     MRNA, 
			     $mRNA->start,
			     $mRNA->end,
			     '.',
			     $mRNA->strand < 0 ? '-' : '+',
			     '.',
			     { 'ID'     => $mrna_name,
			       'Name'   => "$pname.tr",
			       'Parent' => "$pname.gene"}];

		    my %mRNAxref;
		    my @m_xrefs;
		    if( $mRNA->has_tag('db_xref') ) { 
			for my $xref ( $mRNA->get_tag_values('db_xref') ) {
			    my ($xref_src,$xref_id) = split(':',$xref);			    
			    $mRNAxref{$xref_src} = $xref_id;
			    push @m_xrefs, &escape($xref);
			}
		    }

		    $x->[8]->{'Dbxref'} = join(",",@m_xrefs);
		    my @note;
		    if( $mRNA->has_tag('note') ) {
			@note= $mRNA->get_tag_values('note');
		    }
		    if( $mRNA->has_tag('product') ) {
			unshift @note, &escape( $mRNA->get_tag_values('product') );
		    }		    
		    
		    $x->[8]->{'Note'} = join(",",&escape(@note));
		    my @alias;
		    for my $t ( qw(protein_id transcript_id synonym) ) {
			if( $f->has_tag($t) ) {
			    # warn("$t --> ", $f->get_tag_values($t), "\n");
			    push @alias, &escape( $f->get_tag_values($t) );
			}
		    }
		    if( @alias ) {
			$x->[8]->{'Alias'} = join(",", @alias);
		    }

		    my (@exon_a,@cds_a);
		    for my $CDS ( $mRNA->get_SeqFeatures ) {
			if( $CDS->primary_tag eq 'CDS' ) {
			    my (@cdsrefs_a,%cdsxrefs);
			    if(! @cdsrefs_a && $CDS->has_tag('db_xref') ) {
				for my $cxref ( $CDS->get_tag_values('db_xref') ) {
				    my ($cxref_src,$cxref_id) = split(':',$cxref);
				    $cdsxrefs{$cxref_src} = $cxref_id;
				    push @cdsrefs_a, &escape($cxref);
				} 
			    }
			    
			    my $codon_start = 1;
			    if( $CDS->has_tag('codon_start') ) {
				($codon_start) = $CDS->get_tag_values('codon_start');
			    }
			    my $cdsct = 1;
			    for my $cdsloc ( sort { 
				($a->strand * $a->start) <=> 
				    ($b->strand * $b->start)} 
					     $CDS->location->each_Location ) {
				my $id;
				if( exists $cdsxrefs{'GI'} ) {
				    $id = $cdsxrefs{'GI'};
				}
				if( ! defined $id ) {
				    $id = sprintf("%s_cds%05d",$pref,$cdsnum++);
				} else {
				    $id = sprintf("%s.cds%d",$id,$cdsct);
				}
				
				my ($cds_start,$cds_end) = ($cdsloc->start, 
							    $cdsloc->end);
				if( $codon_start > 1 && $cdsct == 1 ) {
				    if( $cdsloc->strand < 0 ) {
					$cds_end -= ($codon_start-1);
				    } else {
					$cds_start += ($codon_start-1);
				    }
				}
				if( $cds_start == $cds_end ) {
				 #    $cds_end+=1;
				}
				push @cds_a, [$seq->display_id,
					      $SRC,
					      CDS, 
					      $cds_start,
					      $cds_end,
					      '.',
					      $cdsloc->strand < 0 ? '-' : '+',
					      '.',
					      { 'ID' => $id,
						'Parent' => $mrna_name}];		
				if( @cdsrefs_a ) {
				    $cds_a[-1]->[8]->{'Dbxref'} = join(",",@cdsrefs_a);
				}
				
				for my $copy ( qw(protein_id codon_start transl_table
						  inference translation) ) {
				    if( ! exists $x->[8]->{$copy} && 
					$CDS->has_tag($copy) ) {
					($x->[8]->{$copy}) = $CDS->get_tag_values($copy);
				    }
				}
				
				$cdsct++;
			    }
			} else {
			    push @exon_a, [$seq->display_id,
					   $SRC,
					   EXON, 
					   $CDS->start,
					   $CDS->end,
					   '.',
					   $CDS->strand < 0 ? '-' : '+',
					   '.',
					   { 'Parent' => $mrna_name}];
			}
		    }
		    push @ALL, $x;
		    push @ALL, @exon_a, @cds_a;
		    $mrnact++;
		}
	    }
	}
    }
    for my $feature ( @ALL ) {
	my $lastcol = pop @$feature;
	my @lastcol_n;
	for my $field ( qw(ID Parent) ) {
	    if( exists $lastcol->{$field} ) {
		push @lastcol_n, "$field=".$lastcol->{$field};
		delete $lastcol->{$field};
	    }
	}	
	for my $k ( keys %{$lastcol}) {
	    push @lastcol_n, sprintf("%s=%s",$k,$lastcol->{$k});
	}
	print $outfh join("\t", @$feature, join(";", @lastcol_n)),"\n";
    }
    #if( $outfh && @ALL && @seqs) {
#	my $out = Bio::SeqIO->new(-format => 'fasta',
#				  -fh     => $outfh);
#	$out->write_seq(@seqs);
#    }
}
$unflattener->report_problems;
sub get_species_gbkfile {
    my ($file) = @_;
    open(SPGREP, "grep '^  ORGANISM' $file | ") || die $!;
    my $species = <SPGREP>;
    close(SPGREP);
    $species =~ s/^\s+ORGANISM\s+//;
    chomp($species);
    
    # -- #
    # Hard code for the current data now

    $species =~ s/var\.\s+//;    
    $species =~ s/(CBS|NRRL)\s(\S+)/$1-$2/;
    
    # -- #
    my $firstl = lc substr($species,0,1);
    substr($species,0,1,$firstl);
    
    $species =~ s/\s+/_/g;
    my ($sp,$gen) = split(/\_/,$species);
    my $pref = substr($sp, 0,1).substr($gen,0,3);
    if( $pref =~ /cneo|cgat|scer/ ) {
	my ($g,$s,@rest) = split(/_/,$species);
	$pref .= "_".pop @rest;
    }
    return ($species,$pref);
}

sub uniqlst { 
    my %x;
    return grep { ! $x{$_}++}  @_;
}

sub escape {
    
    for my $value ( @_) {
	if(  defined $value && length($value) ) { 
	    if ($value =~ /[^a-zA-Z0-9\,\;\=\.:\%\^\*\$\@\!\+\_\?\-]/) {
		$value =~ s/\t/\\t/g;       # substitute tab and newline 
		# characters
		$value =~ s/\n/\\n/g;       # to their UNIX equivalents
		
# Unescaped quotes are not allowed in GFF3
#                   $value = '"' . $value . '"';
	    }
	    $value =~ s/([\t\n\r%&\=;,])/sprintf("%%%X",ord($1))/ge;
	} 
    }
    return @_;
}
