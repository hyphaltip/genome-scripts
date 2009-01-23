#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Long;
my $debug = 0;

GetOptions(
	   'v|verbose!' => \$debug,
	   );
my $genome = shift;

my %codes;
my $code_ct = 0;
my %lengths;

my $in = Bio::SeqIO->new(-format => 'fasta',
			 -file   => $genome);
my %dbh;
while ( my $seq = $in->next_seq ) {
  my $seqid = $seq->display_id;
  if ( $seqid =~ /contig(7\.\d+)/ ) {
    $seqid = "ncra_$1"; # simplify
  }
  $lengths{$seqid} = $seq->length;
}
my $i = 0;
while ( <> ) {
  next if( /^\#/ || /^\s+$/);
  my ($seqid,$src,$type,$start,$end) = split;
  if ( $seqid =~ /contig(7\.\d+)/ ) {
    $seqid = "ncra_$1"; # simplify
  }
  next unless defined $lengths{$seqid};

  unless( exists $codes{$type} ) {
    $codes{$type} = 2**$code_ct++;
  }
  for my $bp ( $start..$end ) {
    if ( ! exists $dbh{$seqid}->{$bp} ||
	 ! ($dbh{$seqid}->{$bp} & $codes{$type}) ) {
      $dbh{$seqid}->{$bp} += $codes{$type};
    }
  }
  warn("$start..$end $type\n") if $debug;
  last if $i++ > 10 && $debug;
}

#if ( $debug) {
  for my $code ( sort { $codes{$a} <=> $codes{$b} } keys %codes ) {
    warn "\t$code $codes{$code}\n";
  }
#}

my %sum;
while ( my ($seqid,$dbh) = each %dbh ) {
  for (my $i =1; $i < $lengths{$seqid}; $i++ ) {
    my $val = 0;
    if ( exists $dbh->{$i} ) {
      $val = $dbh->{$i};
    }
    if ( $val == 0 ) {
      $sum{'NONE'}++;
    } else {
      while ( my ($type,$code) = each %codes ) {
	if ( $val & $code) {
	  $sum{$type}++;
	}
      }
    }
  }
}
print join("\t",qw(TYPE COUNT)),"\n";
for my $type ( sort keys %sum ) {
  print join("\t", $type, $sum{$type}),"\n";
}
