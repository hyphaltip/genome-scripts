use strict;
use warnings;
use Getopt::Long;
my $onlymaker = 1;

GetOptions(
	'onlymaker!' => \$onlymaker,
	);

my $cds = 1;
my $utr3 = 1;
my $utr5 = 1;
while(<>) { 

 if( /^##FASTA/ || /^>/ ) {
   last;
 } elsif( /^\#/) {
   print $_;
   next;
 }
 my @line = split(/\t/,$_);
 if( ! defined $line[1] ) {
  die("line was $_");
 }
 next if $onlymaker && $line[1] ne 'maker';
 if( s/:cds;Parent=/:cds:$cds;Parent=/ ) {
   $cds++; 
 } elsif( s/:three_prime_utr;Parent=/:three_prime_utr:$utr3;Parent=/ ) {
  $utr3++;
 } elsif( s/:five_prime_utr;Parent=/:five_prime_utr:$utr3;Parent=/ ) {
  $utr5++;
 }
 chomp;
 my @row = split(/\t/,$_);
 if( $row[2] eq 'exon' ) {
     my %group = map { split(/=/,$_) } split(/;/,pop @row);
     if( $group{'Parent'} =~ /,/ ) {
	 my @isoforms = split(/,/,$group{'Parent'});
	 my $count = 1;	 
	 my $id = $group{'ID'};
	 for my $p ( @isoforms ) {
	     $group{'ID'} = $id .".".$count++;
	     $group{'Parent'} = $p;	 
	     print join("\t", @row, join(";", 
					 map { sprintf("%s=%s",$_,$group{$_}) } 
					 qw(ID Parent)),';'),"\n";
	 }
     } else {
	 print $_,"\n";
     }
 } else {
     print $_,"\n";   
 }
}
