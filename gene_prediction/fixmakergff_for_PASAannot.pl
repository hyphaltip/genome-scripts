my $cds = 1;
my $utr3 = 1;
my $utr5 = 1;
while(<>) {
 if( s/:cds;Parent=/:cds:$cds;Parent=/ ) {
   $cds++; 
 } elsif( s/:three_prime_utr;Parent=/:three_prime_utr:$utr3;Parent=/ ) {
  $utr3++;
 } elsif( s/:five_prime_utr;Parent=/:five_prime_utr:$utr3;Parent=/ ) {
  $utr5++;
 }
 print;
}
