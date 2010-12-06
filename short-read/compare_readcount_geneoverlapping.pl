#!/usr/bin/perl -w
use strict;
use File::Spec;

# from genome-scripts/seqfeature/find_overlapping_genes.pl
my $overlapping = shift @ARGV || die "need overlapping gene info"; 

# output file(s) from genome-scripts/short-read/SAM_summarize_fwindow.pl
my @read_counts = @ARGV;

open(my $fh => $overlapping) || die "Cannot open $overlapping: $!";
my $i = 0;
my %header = map { $_ => $i++ } split(/\s+/,<$fh>);
my %overlapping_genes;
while(<$fh>) {
    my @row = split;

    push @{$overlapping_genes{ $row[ $header{'MRNA_LEFT'} ] }->{OVERLAP_TYPE} },
    $row[ $header{'OVERLAP_TYPE'}];
    
    push @{$overlapping_genes{ $row[ $header{'MRNA_LEFT'} ] }->{STRAND_MATCH} },
    $row[ $header{'STRAND_MATCH'}];

    $overlapping_genes{ $row[ $header{'MRNA_LEFT'}] }->{MATCH} = $row[ $header{'MRNA_RIGHT'} ];

# -- the overlapping gene also has this information
    push @{$overlapping_genes{ $row[ $header{'MRNA_RIGHT'} ] }->{OVERLAP_TYPE} },
    $row[ $header{'OVERLAP_TYPE'}];

    push @{$overlapping_genes{ $row[ $header{'MRNA_RIGHT'} ] }->{STRAND_MATCH} },
    $row[ $header{'STRAND_MATCH'}];

    $overlapping_genes{ $row[ $header{'MRNA_RIGHT'}] }->{MATCH} = $row[ $header{'MRNA_LEFT'} ];
}
close($fh);

for my $datfile ( @read_counts ) {
    my (undef,undef,$f) = File::Spec->splitpath($datfile);
    next unless ( $datfile =~ /(\S+)\.counts\.dat$/);
    my $stem = $1;
    open($fh => $datfile) || die $!;
    $i = 0;
    my %dat;

    my %d_header = map { $_ => $i++ } split(/\s+/,<$fh>);
    open(my $ofh => ">$stem.overlap_status_readcount_length.dat") || die "cannot open for writing: $!";
    open(my $tofh => ">$stem.overlap_status_readcount.dat") || die "cannot open for writing: $!";
    print $tofh join("\t", qw(FEATURE READCOUNT OVERLAPS OVERLAPS_MRNA OVERLAP_TYPE OVERLAP_STRAND)),"\n";
    
    print $ofh join("\t", qw(FEATURE READLENGTH READCOUNT OVERLAPS OVERLAPS_MRNA 
			     OVERLAP_TYPE OVERLAP_STRAND)),"\n";
    
    while(<$fh>) {
	my @row = split;
	my $total = $row[ $d_header{'SAME'}] + $row[ $d_header{'OPP'}];
	$dat{$row[ $d_header{'FEATURE'} ]} += $total;
	my $has_overlap = exists $overlapping_genes{ $row[ $d_header{'FEATURE'} ] } ? 1 : 0;
	print $ofh join("\t", ( $row[ $d_header{'FEATURE'} ],
				$row[ $d_header{'LENGTH'} ],
				$total,
				$has_overlap,
				$has_overlap ? $overlapping_genes{ $row[ $d_header{'FEATURE'} ] }->{MATCH} : 0,
				$has_overlap ? join(";",@{$overlapping_genes{ $row[ $d_header{'FEATURE'} ] }->{OVERLAP_TYPE}}) : 0,
				$has_overlap ? 
				join(";",@{$overlapping_genes{ $row[ $d_header{'FEATURE'} ] }->{STRAND_MATCH}}) : 0 )),"\n";
    }
    
    for my $gene ( sort { $dat{$b} <=> $dat{$a} } keys %dat ) {
	my $has_overlap = exists $overlapping_genes{$gene} ? 1 : 0;
	print $tofh join("\t", ( $gene, $dat{$gene},
				 $has_overlap,
				 $has_overlap ? $overlapping_genes{$gene}->{MATCH} : 0,
				 $has_overlap ? join(";",@{$overlapping_genes{$gene}->{OVERLAP_TYPE}}) : 0,
				 $has_overlap ? 
				 join(";",@{$overlapping_genes{$gene}->{STRAND_MATCH}}) : 0 )),"\n";
    }
    close($ofh);
    close($tofh);
}
