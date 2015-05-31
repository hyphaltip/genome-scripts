#!env perl
use strict;
use warnings;
use Getopt::Long;

# usage
# generate a lookup table for SRA/ERA/DRP to strain designation
# so we can use strain names in downstream analyses

# parameters - 
# --info a dir containing the SraRunTable.txt downloaded from SRA
# associated with a run 
my $infodir = 'info';

GetOptions('i|info:s' => \$infodir);
#skipping these tags
# Organism_s 
my @cols = qw(SRA_Study_s SRA_Sample_s Sample_Name_s LibraryLayout_s);

opendir(INFO,$infodir) || die "cannot open dir: $infodir: $!";
foreach my $file ( readdir(INFO) ) {
    next unless $file =~ /(\S+)_SraRunTable.txt/;
    # assume files are all end in _SraRunTable.txt
    open(my $fh => "$infodir/$file" ) || die "cannot open $infodir/$file: $!";
    my $header = <$fh>;
    chomp($header);
    my @header = split(/\t/,$header);
    my $i =0;
    my %header2col = map { $_ => $i++ } @header;
    warn sprintf "file: $file cols for @cols\n %s\n",
    join(' ',map { $header2col{$_} } @cols);
    while(<$fh>) {
	chomp;
	my @row = split(/\t/,$_);
	print join("\t",map { $row[ $header2col{$_} ] } @cols), "\n";
	
    }
}
