#!/usr/bin/perl -w

=head1 NAME 

process_exonerate_gff3

=head1 USAGE

process_exonerate_gff3 file1.e2g.EXONERATE > out.gff

=head1 DESCRIPTION

Turns EXONERATE gff output into GFF2 for Apollo use.

I am using 
You need to have run exonerate with at least the following options

 --showtargetgff yes
 --showvulgar yes
 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n"

=head1 AUTHOR

Jason Stajich, jason-at-bioperl-dot-org

=cut

use strict;
use IO::String;
use Getopt::Long;
use Env;
use File::Spec;
my $state = 0; my $ct = 0;
my $buffer = '';
my (%dat,%lengths,%alnlens,%counter,$vulgar,$lastvulgar);
my $outfile;
my $type ='mRNA';
GetOptions(
           'h|help'        => sub { exec('perldoc', $0); exit },
	   'o|out:s'       => \$outfile,
	   );

my $outfh;
if( defined $outfile ) {
    open($outfh, ">$outfile") || die ("cannot open $outfile for writing\n");
} else { 
    $outfh = \*STDOUT;
}

my $first = 1;
while(<>) {
if( $state == 0 && s/^vulgar:\s+// ) {
	chomp($vulgar = $lastvulgar = $_);
	chomp($lastvulgar);
    }
    if( $state > 0 ) {	 
	if( $state == 2 ) {
	    if ( s/^vulgar:\s+//) {
		$lastvulgar = $vulgar;
		chomp($vulgar = $_);
	    }
	    
	    my ($qname,$qstart,$qend, $qstrand,
		$hname,$hstart,$hend,$hstrand,$overall_score,
		@vulgar_process_query) = split(/\s+/,$lastvulgar);
	    my $running_query = $qstart+1; # 1 based system
	    my @cds_triples;
	    while( @vulgar_process_query ) {
		my $label = shift @vulgar_process_query;
		my $qlen  = shift @vulgar_process_query;
		my $hlen  = shift @vulgar_process_query;
		push @cds_triples, [$label,$qlen,$hlen];
	    }
	    my $iostring = IO::String->new($buffer);
	    
	    my ($gene_name,$gene);
	    my $length;
	    while(<$iostring> ) {
		my ($seqid,$srctag,$ptag, $start,$end,$score,$strand,
		    $frame,$groupfield) = split(/\s+/,$_,9);		
		$srctag =~ s/\:/_/g;
		# this is gff2 for exonerate
		$groupfield =~ s/\s+;\s+$//; # drop last field
		my %fields;
		for my $set ( split(/\s+;\s+/, $groupfield) ) {
		    my ($key,$value) = split(/\s+/,$set,2);
		    $fields{$key} = $value;		
		}
		if( $ptag eq 'gene' ) {
		    $length = $lengths{$seqid};		    
		    if(  ! defined $length ) {
			die("unknown length for '$seqid'\n");
		    }
		    $gene = $fields{'sequence'};
		    if ( defined $fields{'gene_orientation'}) {
			if( $dat{$gene} ) {
			    $gene .= $counter{$gene};
			}
			$dat{$gene} = $fields{'gene_orientation'};
			$dat{$gene} = 0 if $dat{$gene} eq '.';
		    }
		    $gene_name = $gene;
		    if( $counter{$gene_name}++ ) {
			$gene_name .= ".$counter{$gene_name}";
		    }
		    next;
		}
		next unless $ptag eq 'exon';
		if( $strand eq '-' ) {
		    my $s = $length - $end + 1;
		    my $e = $length - $start +1;
		    ($start,$end) = ( $s,$e );
		}
		if( $gene && 
		    $dat{$gene} &&
		    $dat{$gene} eq '-' ) 
		{
		    $strand = '-';
		}
		my $match_begin = $running_query;
		my $match_end   = $running_query;
		while( @cds_triples ) {
		    my $triple = shift @cds_triples;
		    last if $triple->[0] eq 'I';
		    $match_end += $triple->[1];
		}
		$running_query = $match_end;
		print $outfh join("\t",
				  $gene_name,$type,'exon',
				  $start,$end,
				  $overall_score, 
				  $strand, $frame, $gene_name, 
				  $match_begin,$match_end,
				  ),"\n";
		
	    }
	    $state = 0;
	    $buffer = '';
	    %dat = ();
	    
	    next;
	} elsif(/^\#/) { next; }
	if( /^>(\S+)\s+length=(\d+)\s+alnlen=(\d+)/ ) {
	    $lengths{$1} = $2;
	    $alnlens{$1} = $3;
	    $ct++;
	} else {	
	    $buffer .= join("\t",split(/\s+/,$_,9));
	}
	if( $ct == 2 ) { 
	    $state = 2;
	    $ct = 0;
	}
    } elsif( /^\# --- START OF GFF DUMP ---/ ) {
	$state = 1; 
	$ct =0;
	$buffer = '';
    }
}
