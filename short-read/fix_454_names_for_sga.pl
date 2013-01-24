#!/usr/bin/perl -w
use strict;
use warnings;

=head1 NAME

fix_454_for_sga.pl

=head1 USAGE 

 perl fix_454_for_sga.pl FILE > NEWFIL

=head1 DESCRIPTION

Fix the read names to be /A /B instead of ending in a,b

=cut1

while(<>) {
 my $id = $_;
 if( $id =~ /^(\@\S+)([ABab])(\s+.*)/ ) {
  $id = sprintf("%s/%s%s",$1,uc($2),$3);
  if( $id !~ /\n$/ ) { $id .= "\n"; }
 } else {
  chomp($id);
  die("cannot process ID - $id doesn't match any pattern\n");
 }
 my $seq = <>;
 my $dsc = <>;
 my $qual = <>;
 print $id, $seq, "+\n", $qual;
}
