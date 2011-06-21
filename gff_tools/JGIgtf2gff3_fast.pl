#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;

use constant CDSEXON => 'cds';
use constant EXON => 'exon';
use constant MRNA => 'mRNA';
use constant GENE => 'gene';
my $verbose = 0;

GetOptions('v|verbose!'     => \$verbose,
	   );
my $src_string = "JGI";

for my $file ( @ARGV ) {
    my $stem = $file;
    $stem =~ s/\.g[ft]f//;
    open(my $fh=> $file) || die $!;
    my (%genes);
    my %lookup;
    while(<$fh>) {
	my ($seqid,$src,$type,$start,$end,$score,
	    $strand,$frame,$lastcol) = split(/\t/,$_,9);
	next unless $type eq 'CDS' || $type eq 'exon';
	
	if( ! defined $seqids{$seqid} ) {
	    warn("no sequence for $seqid\n");
	    exit;
	} elsif ( $seqids{$seqid} < $end ) {
	    warn("seq length of $seqid is $seqids{$seqid} is less than annotation end ($end)\n");
	    exit;
	}
	
	my ($geneid,$transcriptid);
	if( $lastcol =~ /(?:name|gene_id)\s+\"([^\"]+)\";/ ) {
	    $geneid = $1;
	}
	
	if( $lastcol =~ /featureId\s+(\d+)/ ) {
	    $transcriptid = "$prefix\_t_$1";
	}
	if( $lastcol =~ /transcriptId\s+(\d+)/ ) {
	    $transcriptid = "$prefix\_t_$1";
	} elsif( $lastcol =~ /transcript_id\s+\"([^\"]+)\"/ ) { #"
	    $transcriptid = "$prefix\_t_$1";
        }
	
	if( defined $transcriptid && ! defined $lookup{$geneid}) {
	    $lookup{$geneid} = $transcriptid;
	}
	if( ! defined $transcriptid ) {
	    $transcriptid = $lookup{$geneid};
	}
	if( ! defined $geneid || ! defined $transcriptid ) {
	    warn("cannot find geneid or transcriptid in $lastcol\n");
	    next FILE;
	}
	my $typetag = $type eq 'CDS' ? CDSEXON : EXON;

	push @{$genes{$geneid}->{$transcriptid}}, 
	[$seqid,$typetag, $start,$end,$score,$strand,$frame];
    }
    open(my $out => ">$$stem.gff3") || die $!;
    print $out "##gff-version 3\n";
    print $out "#date-prepared ".localtime(time)."\n";
    my %num= ('exon' => 1, 'cds' => 1);
    for my $gene ( keys %genes ) {
	my ($gmin, $gmax,$gstrand,$seqid_all);
	my @lines;
	for my $transcript ( keys %{$genes{$gene}} ) {
	    my ($tmin,$tmax,$tstrand);
	    for my $CDS ( @{$genes{$gene}->{$transcript}} ) {
		my ($seqid,$type,$start,$end,$score,$strand,$frame) = @$CDS;
		($start,$end) = sort { $a <=> $b } ($start,$end);
		$gmin = $start if ! defined $gmin || $gmin > $start;
		$gmax = $end if ! defined $gmax || $gmax < $end;
		$gstrand = $strand;

		$tmin = $start if ! defined $tmin || $tmin > $start;
		$tmax = $end if ! defined $tmax || $tmax < $end;
		$tstrand = $strand;
		
		push @lines, join("\t", $seqid, $src_string, $type,
				  $start, $end, $score, $strand, $frame,
				  join(';',
				       sprintf("ID=%s_%s%05d",$prefix,lc($type),
					       $num{$type}++),
				       sprintf('Parent=%s',$transcript)));
		$seqid_all = $seqid;
	    }
	    unshift @lines, join("\t", $seqid_all, $src_string, MRNA,
				 $tmin, $tmax, '.', $tstrand, '.',
				 join(';',
				      sprintf("ID=%s;Name=%s",$transcript,$transcript),
				      sprintf('Parent=%s',$gene)));
	}
	print $out join("\t", $seqid_all, $src_string, GENE,
			$gmin, $gmax, '.', $gstrand, '.',
			join(';',
			     sprintf("ID=%s",$gene),
			     sprintf("Name=%s",$gene))
			), "\n", join("\n", @lines), "\n";
    }
}

