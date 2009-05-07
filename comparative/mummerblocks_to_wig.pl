#!/usr/bin/perl -w
use strict;

=head1 NAME

=head1 USAGE

=head1 AUTHOR

=cut
# take blocks from mummercoords_to_blocks and 
# their aligned files and convert into WIG file

use Bio::SeqIO;
use Bio::SearchIO;
use Bio::AlignIO;
use Getopt::Long;

my ($db,$dir,$out);
my ($afile,$aformat,$seqformat) = qw(align.GGSEARCH fasta fasta);
my $debug = 0;
my $window = 100;

my $name = 'MUMMER';
my $desc = 'MUMMER Percent Identity';
GetOptions('v|verbose!' => \$debug,
	   'f|format:s' => \$aformat,
	   'a|align:s'  => \$afile,
	   'db:s'       => \$db,
	   'd|dir:s'    => \$dir,
	   'o|out:s'    => \$out,
	   'w|window:i' => \$window,
	   'desc:s'     => \$desc,
	   'n|name:s'   => \$name,
	   );

unless( defined $dir ) {
    die "Need a directory of blocks with -d/--dir\n";
}
unless( defined $db ) {
    die "Need a genome database file (FASTA) with -db\n";
}
unless( defined $out ) {
    $out = "$dir.$db.wig"
}
my $refidx = 1; # reference sequence is the db in the GGSEARCH queries
my %genome;
my %lengths;
my $in = Bio::SeqIO->new(-format => $seqformat,
			 -file   => $db);
while( my $seq = $in->next_seq ) {
    my $id        = $seq->display_id;
    $lengths{$id} = $seq->length;
    $genome{$id}->[int($seq->length / $window)] = 0;    
}
opendir(DIR, $dir) || die "cannot open $dir: $!";
for my $subdir ( readdir(DIR) ) {
    next unless ( $subdir =~ /^(\d+)$/ );
    warn("block $subdir\n") if $debug;
    my $alnfile = File::Spec->catfile($dir,$subdir,$afile);
    unless( -f $alnfile ) {
	warn("Skipping $subdir, no $afile\n");
	next;
    } else {
	my $in = Bio::SearchIO->new(-format => $aformat,
				   -file   => $alnfile);
	
	R: while( my $r = $in->next_result ) {
	    my $qname = $r->query_name;
	  
	    my ($qid,$qstart,$qend);
	    if( $qname =~ /(\S+)_(\d+)\-(\d+)$/ ) {
		($qid,$qstart,$qend) = ($1,$2);
	    } else {
		warn("cannot parse $qname for location\n");
		next;
	    }
	    if( my $h = $r->next_hit ) {
		my $hname = $h->name;
		my ($hid,$hstart,$hend);
		if( $hname =~ /(\S+)_(\d+)\-(\d+)$/ ) {
		    ($hid,$hstart,$hend) = ($1,$2);
		} else {
		    warn("cannot parse $hname for location\n");
		    next R;
		}
		if( my $hsp = $h->next_hsp ){
		    my $aln = $hsp->get_aln;
		    my ($window_start) = int($hstart / $window);
		    my $offset = ($hstart % $window);
		    warn("$hid: $hstart $offset $window_start \n") if $debug;
		    my @ids;
		    my @seqs = map { push @ids, $_->id;
				     $_->seq } $aln->each_seq;
		    my $c = length($seqs[$refidx]);
		    # remove the gaps from the ref since we want ref coords 
		    while( ($c = rindex($seqs[$refidx],'-',$c)) > 0) {
			substr($seqs[$refidx],$c,1,'');
			substr($seqs[1-$refidx],$c,1,'');
		    }
		    if( $offset > 1 ) {
			my @preblock = map { substr($_,0,$offset) } @seqs;
			warn("preblock is \n",join("\n",@preblock),"\n") 
			    if $debug;
			$genome{$hid}->[$window_start++] = &percent_id(@preblock);
		    }
		    my $len = $aln->length;
		    for( my $i = $offset; $i < $len; $i+= $window ) {
			my $window_cut = $window;
			
			if( $i+$window >= $len ) {
			    $window_cut = $len - ($i+$window);
			}
			my @slice = map { substr($_,$i,$window) } @seqs;
			warn("i is $i slice is \n",join("\n",@slice),"\n") 
			    if $debug;
			$genome{$hid}->[$window_start++] = &percent_id(@slice);
		    }	    
		} else {
		    warn("No HSP for $qname/$hname\n");
		}
	    } else {
		warn("No hit for $qname\n");
	    }
	}
    }
    last if $debug;
}

open(my $fh => ">$out") || die "$out; $!";
printf $fh "track type=wiggle_0 name=\"%s\" description=\"%s\"\n",
    $name, $desc;

for my $chrom ( sort { $lengths{$b} <=> $lengths{$a} } 
		keys %genome ) {
    printf $fh "fixedStep chrom=%s start=%d step=%d\n",
    $chrom, 1, $window;
    
    for my $w ( @{$genome{$chrom}} ) {
	print $fh $w;
    }
}
