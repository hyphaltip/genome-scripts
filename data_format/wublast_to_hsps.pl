#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my (%group,$groupcounter,%fields);
my ($src,$type,$score) = qw(WUBLAST match pcident);
my ($minid,$minlen) = (0,0);
GetOptions(
	   's|src:s' => \$src,
	   't|type:s'=> \$type,
	   'score:s' => \$score, # name of scoring field
	   'm|minid:i' => \$minid,
	   'l|minlen:i'=> \$minlen,
	   );

my @rows;
while(<>) {

    if( /^\#/ ) {
	# process some header info
	if( s/\# Fields:\s+//) {
	    my $i = 0;
	    %fields = map { $_ => $i++ } split;
	} 
	next;
    }    
    chomp; # drop trailing \n since we spliting by tab
    my @line = split(/\t/,$_);
    my $groupid;
    my $target = $line[ $fields{'sid'} ]; # target name

    # now build up an ID name based on either the unique links groups
    # or just assign everything from the same target as the group name
    if( exists $fields{'links'} ) {
	my $link = $line[ $fields{'links'} ];
	$link =~ s/[\(\)]//g;
	unless( exists $group{$target .":".$link} ) {
	    $group{$target .":".$link} = $groupcounter++;
	}
	$groupid = $group{$target .":".$link};	
    } else {
	$groupid = $group{$target}++ ? $groupcounter : $groupcounter++;
    }
    my ($start,$end,$strand) = ($line[ $fields{'qstart'}], 
				$line[ $fields{'qend'}],'+');
    if( $start > $end ) { # flip-flop so at start < end and strand is flipped
	($start,$end,$strand) = ($end,$start,'-');
    }
    next if ( $line[ $fields{$score}] < $minid ||
	      $line[ $fields{'alignlen'}] < $minlen );
    push @rows, [$line[0],
		 $src, $type,	       
		$start,$end,
		 $line[ $fields{$score} ],	       
		 $strand,
		 '.', # could be phase from translated search?
		 sprintf("ID=Match%s;Name=%s.%d;Target=%s %d %d",
			 $groupid, $target, $groupid,$target, 
			 $line[ $fields{'sstart'}],
			 $line[ $fields{'send'}])];    
}
for my $row ( sort { $a->[-1] cmp $b->[-1] } @rows ) {
    print join("\t", @$row), "\n";
}
