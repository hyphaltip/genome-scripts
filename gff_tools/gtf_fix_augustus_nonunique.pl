while(<>) {
    my @line = split(/\t/,$_);
    my %dat;
    
    for my $field ( qw(gene_id transcript_id) ) {	
	if( $line[-1] =~ s/$field\s+\"(\S+)\"/$field "$line[0]-aug.$1"/ ) {
	}
    }
  print join("\t", @line),"\n";
}
