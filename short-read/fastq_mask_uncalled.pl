
while(<>) {
    my $head = $_;
    my $seq  = <>;
    my $desc = <>;
    my $qual = <>;
    if( substr($head,0,1) ne '@' ||
	substr($desc,0,1) ne '+' ) {
	die("register for parsing is off");
    }
    $seq =~ tr/\.\-/NN/;
    print $head, $seq, "+\n",$qual;
}
