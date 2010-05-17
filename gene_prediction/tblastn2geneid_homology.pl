#!/usr/bin/perl -w

# tblastn -links -group
use strict;
use Getopt::Long;
use File::Spec;

my %compress = ('gz'   => 'zcat',
		'bz2'  => 'bzcat',
		);

my $informat = 'wutab'; # wublast mformat 3 is default

my $gffver = 3;
my $output;
my $debug = 0;
my $best_only = 1;
my $evalue = '1e-25';
GetOptions(
	   'o|out:s'         => \$output,
	   'f|format:s'      => \$informat,
	   'v|debug!'        => \$debug,
	   'best!'           => \$best_only,
	   'e|evalue:s'      => \$evalue,
    );


my $out;
if($output) {
    open($out,">$output") || die $!;
} else {
    $out = \*STDOUT;
}

for my $file ( @ARGV) {
    my $fh;
    if( $file =~ /\.(gz|bz2)$/ ) {
	open( $fh => "$compress{$1} $file |") || die $!;
    } else {
	open( $fh => $file) || die $!;
    }
    if( lc($informat) eq 'wutab') {
	my %lastseen;
	my @hsps;
	my $i = 0;

	while(<$fh>) {
	    chomp;
	    next if(/^\#/);
	    my $r = $_;
	    my ($q,$h,$e,$N,$Sprime,$S,
		$alignlen,$nident,$npos,
		$nmism, $pcident,$pcpos,$qgaps,$qgaplen,
		$sgaps,$sgaplen,
		$qframe,$qstart,$qend,
		$sframe,$sstart,$send,
		$group, $links) = split(/\t/,$_);

#	    next if $evalue && $e > $evalue;
#	    warn("group is $group links are $links for $q $h $e\n") if $debug;
	    $links =~ s/[\(\)]//g;

	    if( @hsps &&
		( $lastseen{'query'} ne $q ||
		  $lastseen{'hit'}   ne $h ) ) {
		for my $grp ( keys %{$lastseen{'groups'}} ) {
		    #warn("group is ",join(",",split(/\-/,$group))," $group\n");
		    #warn("$q $h $sstart, $send ", $#hsps," $grp\n");
		    &make_pair($lastseen{'query'},
			       $lastseen{'hit'},
			       [ map { $hsps[$_ - 1] } split(/\-/,$grp)]);
#		    last if $best_only;
		}
		@hsps = ();
		$lastseen{'groups'} = {}; # reset the groups    
		#last if $debug && $i++ > 3;
	    }
	    $lastseen{'groups'}->{$links}++;
	    ($sstart,$send) = sort { $a <=> $b} ($sstart,$send);
	    push @hsps, { qstart  => $qstart,
			  qend    => $qend,
			  hstrand => substr($sframe,0,1),
			  hstart  => $sstart,
			  hend    => $send,
			  pid     => $pcident,
		      };
	    $lastseen{'query'} = $q;
	    $lastseen{'hit'}   = $h;
	}
	if( @hsps ) {
	    for my $grp ( keys %{$lastseen{'groups'}} ) {
		&make_pair($lastseen{'query'},
			   $lastseen{'hit'},
			   [ map { $hsps[$_ - 1] } split(/\-/,$grp)]);
	    }
	}	
    }
}


sub make_pair {
    my ($q,$h,$hsps) = @_;
    my $qname = $q;

    for my $hsp ( sort { $a->{hstart} <=> $b->{hstart} } @$hsps) {
	if( ! defined $hsp->{hstart} ) {
	    warn("no hstart for $q $h ",keys %$hsp,"\n");
	    next;
	}
	my $strand = $hsp->{hstrand};
	print $out join("\t", $h, 'TBLASTN', 'HSP', 
			$hsp->{'hstart'},
			$hsp->{'hend'}, 
			$hsp->{'pid'},
			$strand, 
			'.',
			"Target=$q"),"\n";
    }
}
