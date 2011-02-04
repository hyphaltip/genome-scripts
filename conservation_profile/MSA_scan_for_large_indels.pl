#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use Getopt::Long;
my $Min_gapsize = 50;
my $Max_gapsize = 5000;
my $alndir = 'alignments';
my $genome;
my $strand = '+';
my $agp_ref;
my $exe = "sliceAlignment %s %s %s %d %d %s |";
my $debug = 0;
GetOptions('r|ref=s' => \$genome,
	   'agp:s'   => \$agp_ref,
	   'd|dir:s' => \$alndir,
	   'v|verbose!' => \$debug,
	   );
die"must have genome with -r or --ref\n" unless $genome;
$agp_ref ||= "$genome.agp";

# read in the AGP file from Mercator to get the start/stop of chromosomes
open(my $agp => $agp_ref ) || die "$agp_ref: $!";
my %refchroms;
while(<$agp>) {
    my ($chrom,$start,$end) = split;
    $refchroms{$chrom} = [$start => $end];
}

for my $chrom ( sort keys %refchroms ) { 
    my ($start,$end)= @{$refchroms{$chrom}};
    warn("$chrom $start .. $end\n");
 # run sliceAlignment to get the projected MSA
    open( my $fh => sprintf($exe, $alndir,$genome,$chrom,
			    $start,$end,$strand)) || 
				die "$!";    

    my $alnio = Bio::AlignIO->new(-fh => $fh,
				  -format   => 'fasta');
   # parse alignment
    while( my $aln = $alnio->next_aln ) {
	warn("got aln, len ", $aln->length,"\n");
	my %snv;
	my @s = $aln->each_seq;
	my @d =  map { $_->seq } @s;
	# extract the chromosome:start-stopSTRAND from the description field
	my @desc = map { my ($chm,$loc) = split(/:/,$_->description);
			 my ($s,$e,$strd) = ($loc =~ /(\d+)-(\d+)([+-])/);
			 [($chm,$s,$e,$strd)];
		     } @s;
	my $len = $aln->length;
	# walk down the length of the alignment
	for( my $i = 0; $i < $len; $i++ ) {
	    my %alleles;	    
	    for my $s ( @d ) {
	        # save time by just getting 1st character, deleting it, but saving it, there are no seeks long string this way
		my $ch = substr($s,0,1,'');
		$alleles{$ch}++;
	    }
	    next if keys %alleles == 1; # ignore where monomorphic
	    if( exists $alleles{'-'} ) {
		$snv{$i} = 'GAP'; # code GAPs
	    } else { 
		$snv{$i} = join(",",keys %alleles); # store the SNP alleles
	    }
	}
	my @gaps = sort { $a <=> $b } grep { $snv{$_} eq 'GAP' } keys %snv;
	my @snps = sort { $a <=> $b } grep { $snv{$_} ne 'GAP' } keys %snv;
	my @collapse_gaps = &collapse_nums(@gaps); # run the collapse algorithm to get runs of GAPs
	
	print "there are ", scalar @snps, " SNPs and ", scalar @gaps, " total gaps and ", scalar @collapse_gaps, " gap openings\n";	

	for my $indel ( @collapse_gaps ) {
	    # capture the indels, check if they are right size, and project the location from 
	    if( $indel =~ /(\d+)-(\d+)/ ) { # only those non singleton indels
		my ($from,$to) = ($1,$2);
		my ($size_gap) = abs($to - $from);
		if( $size_gap > $Min_gapsize && $size_gap < $Max_gapsize ) {
		   my ($locfrom) = $s[0]->location_from_column($from); # project from aln space to seq coords
		   my ($locto)  =  $s[0]->location_from_column($to); 
		# Report the info about the gap
		    printf "%s aln=%s %d..%d %d\n",
		    $desc[0]->[0],
		    $indel,		    
		    $desc[0]->[1] + $locfrom->start,
		    $desc[0]->[1] + $locto->start,
		    $size_gap,
		}
	    }
	    
	}
	last if $debug;
    }
    last if $debug;
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
