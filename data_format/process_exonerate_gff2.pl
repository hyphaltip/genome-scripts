#!/usr/bin/perl -w

=head1 NAME 

process_exonerate_gff3 -t EST file1.e2g.EXONERATE file2.e2g.EXONERATE > out.gff2

=head1 USAGE

Command line arguments
 -t|type EST or Protein typically, optional string to prefix the feature's group name by
     i.e.  EST:af123456.r1 or Protein:YAL001C
 -g|groupname The name of the grouping field - 'GenePrediction' is used for now
 --showalign - show alignments features HSP and match (not completely finished for 
               reverse strand features)
 
=head1 DESCRIPTION

Turns EXONERATE gff output into GFF2 for Gbrowse use.

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
my $grouping_name = 'GenePrediction';
my (%dat,%lengths,%alnlens,%counter,$vulgar,$lastvulgar);
my $outfile;
my $type = ''; # or EST
my $showalign = 0;
GetOptions(
	   't|type:s'      => \$type,
           'h|help'        => sub { exec('perldoc', $0); exit },
	   'o|out:s'       => \$outfile,
	   'showalign!'    => \$showalign,
	   'g|groupname:s' => \$grouping_name,
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
		    if( defined $fields{$key} ) {
			push @{$fields{$key}}, $value;
		    } else {
			if( $key eq 'Align' ) {
			    $fields{$key} = [$value];
			} else {
			    $fields{$key} = $value;
			}
		    }
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
 		    for my $t ( qw(gene_id sequence) ) {
			delete $fields{$t};
		    }
		    $gene_name = $gene;
		    if( $counter{$gene_name}++ ) {
			$gene_name .= ".$counter{$gene_name}";
		    }
		    $fields{'vulgar'} = $lastvulgar;
		    my (undef,$min,$max) = split(/\s+/,$lastvulgar,4);
		    ($min,$max) = sort { $a <=> $b } ($min,$max);
		    my $genelen = $lengths{$gene};
		    $fields{'Length'}   = $genelen; 
		    $fields{'Coverage'} = int(100 * ($max-$min)/$genelen);
		    
		    if( $srctag =~ /est2genome/ ) {
			$ptag = 'match';
		    } elsif( $srctag =~ /protein2genome/ ) {
			$ptag = 'mRNA';
		    } else {
			warn("unknown srctag '$srctag'\n");
		    }
			
		} elsif( $ptag eq 'similarity') {
		    next unless $showalign;
		    if( defined $fields{'Align'} ) { 
			my @hsps = map { split(/\s+/,$_) } @{$fields{'Align'}};
			$ptag = 'match';
			
			my ($min,$max);
			while ( @hsps ) {
			    my ($hstart, $qstart,$hlen) = splice(@hsps,0,3);
			    my $hlen_mod = $hlen;
			    if( $srctag =~ /protein2genome/i ) {
				$hlen_mod /= 3;
			    }
			    if( $strand eq '+' ) {
				print $outfh join("\t",
						  $seqid,
						  $srctag,
						  'HSP',
						  $hstart,
						  $hstart + $hlen,
						  '.',
						  $dat{$gene},
						  sprintf("Target %s %d %d",
							  $gene_name,
							  $qstart,
							  $qstart+$hlen_mod),
						  ), "\n";
			    } else {
				if( $first ) {
				    warn("reverse strand alignment writing is screwy...\n");
				    $first = 0;
				}
			    }
			    $min = $qstart unless defined $min;
			    $max = $qstart + $hlen_mod;
			}
			$fields{'Target'} = sprintf("$type:%s %d %d",
						    $gene_name, $min, $max,
						    );
			my $genelen = $lengths{$gene};
			$fields{'Length'}   = $genelen; 
			$fields{'Coverage'} = int(100 * ($max-$min)/$genelen);
			delete $fields{'Align'};
		    }
		} elsif( $ptag eq 'splice5' ) {
		    $ptag = 'splice_donor';
		} elsif( $ptag eq 'splice3' ) {
		    $ptag = 'splice_acceptor';
		} elsif( $ptag eq 'frameshift' ) {
		    $ptag ='frameshift';
		} elsif( $ptag eq 'exon' ) {
		    if( $srctag =~ /est2genome/ ) {
			$ptag = 'HSP';
		    } elsif( $srctag =~ /protein2genome/ ) {
			$ptag = 'CDS';
		    } else {
			warn("unknown srctag '$srctag'\n");
		    }
		}
		for my $tag ( qw(intron_id deletions insertions) ) {
		    delete $fields{$tag};
		}
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

		if( length($type) ) {
		    $fields{$grouping_name} = "$type:$gene_name";
		} else {
		    $fields{$grouping_name} = $gene_name;
		}
		my @fields = grep { ! /^ID|Parent|GenePrediction|\Q$grouping_name\E|Target$/ } sort keys %fields;
		for my $master ( $grouping_name, qw(ID Parent Target) ) {
		    if( defined $fields{$master} ) {
			unshift @fields, $master;
			last;
		    }
		}
		print $outfh join("\t",
				  $seqid,$srctag,$ptag,$start,$end,$score, 
				  $strand, $frame,
				  join(" ; ", map { $_ . " ". $fields{$_} } @fields)), "\n";
		
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
