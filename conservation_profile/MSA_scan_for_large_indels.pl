#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Getopt::Long;

my $Min_gapsize = 40;
my $Max_gapsize = 10000;
my $Max_Npercent = 0.50;

my $alndir = 'alignments';
my $ref_genome;
my $mapfile = 'map';
my $namesfile = 'genomes';
my $exe = "sliceAlignment %s %s %s %d %d %s |";
my $sdbExtract = "sdbExport %s.sdb %s %d %d %s |";
my $debug = 0;
GetOptions('r|ref=s' => \$ref_genome,
	   'mapfile:s' => \$mapfile,
	   'n|names:s' => \$namesfile,
	   'd|dir:s' => \$alndir,
	   'v|verbose!' => \$debug,
	   );
die"must have refgenome with -r or --ref\n" unless $ref_genome;

# read in genomes order
open(my $gnfh => "$alndir/$namesfile" ) || die "$alndir/$namesfile: $!";
my $x = 0;
my %genomeorder = map { $_ => $x++ } split(/\s+/,<$gnfh>);
close($gnfh);

die("refgenome $ref_genome unrecognized\n") unless exists $genomeorder{$ref_genome};
# read in the map file from Mercator to get the start/stop of chromosomes
open(my $map => "$alndir/$mapfile" ) || die "$alndir/$mapfile: $!";
my @chunks;
while(<$map>) {
    my ($chunk,@dat) = split;
    # the columns we want are based on the offset (multiples of 3) 
    # and the genome order, eg the ref is the 1st one usually so it 
    # will be offset 1
    my ($chrom,$start,$end,$strand) = splice(@dat, $genomeorder{$ref_genome} * 4,4); 
    next if $start eq 'NA' || $end eq 'NA';
    push @chunks, [ $chrom,$start,$end,$strand];
}

my $gapcounter = 0;
open( my $ofh => ">$ref_genome.gaps.dat") || die $!;
open(my $gapseqs => ">$ref_genome.gap_seqs.fa") || die $!;

print $ofh join("\t", qw(GAPID CHROM ALNPOS SEQPOS LENGTH)),"\n";
for my $chunk ( @chunks ) { 
    my ($chrom,$chr_start,$chr_end,$strand)= @$chunk;
    # run sliceAlignment to get the projected MSA
    my $cmd = sprintf($exe, $alndir,$ref_genome,$chrom,
		      $chr_start,$chr_end,$strand);
    warn("cmd is $cmd\n") if $debug;
    open( my $fh => $cmd ) || die "$cmd: $!";    

    my $alnio = Bio::AlignIO->new(-fh => $fh,
				  -alphabet => 'dna',
				  -format   => 'fasta');
    # parse alignment
    if( my $aln = $alnio->next_aln ) {
	my %snv;
	my @s = $aln->each_seq;

	# figure out which seq is the 'ref' seq that we've requested based on genome name
	# keep the order number so when we search for it later we know which one to remap to
	my $order = 0;
	my $n = 0;
	for my $id ( map { $_->id } @s ) {
	    if( $id eq $ref_genome ) {
		$order = $n;
		last;
	    }
	    $n++;
	}
	my @d =  map { $_->seq } @s;
	# extract the chromosome:start-stopSTRAND from the description field
	my @desc = map { 
	    my $desc = $_->description || '';
	    my ($chm,$loc) = split(/:/,$desc);
	    my ($s,$e,$strd) = (0,0,'');
	    if(defined $loc && $loc =~ /(\d+)-(\d+)([+-])/) {
		($s,$e,$strd) =($1,$2,$3);
	    }
	    [($chm,$s,$e,$strd)];
	} @s;
	my $len = $aln->length;
	# walk down the length of the alignment
	my %allgaps;	    
	for( my $i = 0; $i < $len; $i++ ) {
	    my %alleles;
	    my $si = 0;
	    for my $s ( @d ) {
	        # save time by just getting 1st character, deleting it, but saving it, there are no seeks long string this way
		my $ch = substr($s,0,1,'');
		$alleles{$ch}++;
		if( $si != $order ) {
		    $allgaps{$i}++ if $ch eq '-';
		}
		$si++;
	    }
	    next if keys %alleles == 1;	# ignore where monomorphic
	    if( ! exists $alleles{'-'} ) {
		$snv{$i} = join(",",keys %alleles); # store the SNP alleles
	    }
	}
	my @gaps = sort { $a <=> $b } grep { $allgaps{$_} > 1 } keys %allgaps;
	my @snps = sort { $a <=> $b } keys %snv;
	my @collapse_gaps = &collapse_nums(@gaps); # run the collapse algorithm to get runs of GAPs
	warn("collapse_gaps = @collapse_gaps\n") if $debug;
	my $singleton = 0;
	my @kept_gaps;
	for my $indel ( @collapse_gaps ) {
	    # capture the indels, check if they are right size, and project the location from 
	    if( $indel =~ /(\d+)-(\d+)/ ) { # only those non singleton indels
		my ($from,$to) = ($1,$2);
		my ($size_gap) = abs($to - $from);
		if( $size_gap > $Min_gapsize && $size_gap < $Max_gapsize ) {
		    # project from aln space to original seq coords
		    my ($locfrom) = $s[$order]->location_from_column($from);	
		    my ($locto)  =  $s[$order]->location_from_column($to); 
		    # Report the info about the gap
		    # COLUMNs in the desc are ('chrom','start','end','strand')
		    my $seqcoord_start = $desc[$order]->[1] + $locfrom->start;
		    my $seqcoord_end   = $desc[$order]->[1] + $locto->start;
		    
		    open( my $extract => sprintf($sdbExtract,$ref_genome,$chrom,
						 $seqcoord_start,
						 $seqcoord_end,
						 $strand)) || die "sdbextract: $!";
		    my $indelseq = '';
		    my $record;
		    while(<$extract>) {
			if(/^>(\S+)/) {
			    $record = sprintf(">%d %s:%d..%d\n", $gapcounter++,$chrom,
					      $seqcoord_start, $seqcoord_end);
			} else {
			    $record .= $_;
			    chomp($indelseq .= $_);
			}	    
		    }
		    my $Ncount = ($indelseq =~ tr/Nn/Nn/);
		    next if( ($Ncount / $size_gap) > $Max_Npercent);
		    
		    # print the gap seq to a file
		    print $gapseqs $record;
		    
		    # print the gap info to the report
		    printf $ofh "%d\t%s\t%s\t%d..%d\t%d\n",
		    $gapcounter,
		    $chrom,
		    $indel,
		    $seqcoord_start, $seqcoord_end,
		    $size_gap;

		    push @kept_gaps, $indel;		    
		}
	    } else { 
		$singleton++;
	    }	    
	}
	print "$chrom $chr_start-$chr_end there are ", scalar @snps, " SNPs and ", 
	scalar @gaps, " total gaps and ", scalar @kept_gaps, " gap openings and ",
	$singleton, " singleton indels\n";	

#	last if $debug;
    }    
}

sub collapse_nums {
# This is probably not the slickest connectivity algorithm, but will do for now.
    my @a = @_;
    my ($from, $to, $i, @ca, $consec);
    
    $consec = 0;
    for($i=0; $i < @a; $i++) {
        not $from and do{ $from = $a[$i]; next; };
    # pass repeated positions (gap inserts)
    next if $a[$i] == $a[$i-1];
        if($a[$i] == $a[$i-1]+1) {
            $to = $a[$i];
            $consec++;
        } else {
            if($consec == 1) { $from .= ",$to"; }
            else { $from .= $consec>1 ? "\-$to" : ""; }
            push @ca, split(',', $from);
            $from =  $a[$i];
            $consec = 0;
            $to = undef;
        }
    }
    if(defined $to) {
        if($consec == 1) { $from .= ",$to"; }
        else { $from .= $consec>1 ? "\-$to" : ""; }
    }
    push @ca, split(',', $from) if $from;

    @ca;
}
