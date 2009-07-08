#!/usr/bin/perl -w
print "##gff-version 3\n";
print "##date ".localtime(time),"\n";
my $type = shift || die "need a type\n";
while(<>){
    my ($chrom,$start,$end,$id,$score) = split;
    print join("\t", $chrom,qw(FSEQ), "$type\_peak",$start,$end,$score,'+','.',
	       sprintf("ID=$typePeak_%s",$id)),"\n";
}
