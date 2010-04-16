#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my $dir =  ".";

GetOptions(
	   'd|dir:s' => \$dir,
	   );

my $pattern = qr/(\S+)_(\d+)-(\d+)\.RNAz/;
my $exp_status = 'RNA';
print "##gff-version 3\n";
opendir(DIR, $dir) || die $!;
my $rnazid = 0;
my $ref;
for my $file ( readdir(DIR) ) {

    if( $file =~ $pattern ) {
	my ($chrom,$start,$stop) = ($1,$2,$3);
	# expect the reference name to be encoded in the first 
	# NAME_CHROM info in the file name
	# In coprinus this means the .RNAz files are called
	# Ccin_Chr_10_1009713-1009833.RNAz
	# and the names of sequences in the file are
	# >Ccin
	# >Lbic
	# >Scom 
	# etc
	my ($ref) = split(/_/,$chrom);
#	warn("$file $chrom:$start..$stop\n");
	open(my $fh => "$dir/$file") || die $!;
	my (%dat,$seen_hash,%seqs,%struct);
	while(<$fh>) {
	    if( /\s*(.+):\s+(\S+)/ ) { # matcing the NAME: DATA for RNAz output
		$dat{$1} = $2;		
#		warn("$1 --> $2\n");
	    } elsif( /^\#+\s+RNAz\s+(\S+)\s+\#+/ ) {
		# skip
	    } elsif( /^\#+\s+$/ ) {
		$seen_hash = 1;
	    } elsif( $seen_hash ) {
		if( /^>(\S+)/ ) {
		    my $seqid = $1;
		    $seqs{$seqid} = <$fh>;
		    my $tmp;
		    ($struct{$seqid},$tmp) = split(/\s+/,<$fh>,2);
		    if( $tmp =~ /\(\s*(-?\d+\.\d+)/) {			
			($dat{'MFE'}->{$seqid}) = $1;
		    } 
		}
	    }
	}
	close($fh);
	$dat{'Reading direction'} = 
	    $dat{'Reading direction'} eq 'forward' ? '+' : '-';
	if( $exp_status eq $dat{'Prediction'} ) {
	    print join("\t",
		       $chrom,'RNAz','ncRNA',$start,$stop,
		       $dat{'SVM RNA-class probability'},
		       $dat{'Reading direction'},'.',
		       sprintf("ID=ncRNA%04d;Name=ncRNA%04d;Consensus_MFE=%s;RefMFE=%s;RefStruct=%s",
			       $rnazid,$rnazid,
			       $dat{'Consensus MFE'},
			       $dat{'MFE'}->{$ref},
			       $struct{$ref},
			       )),"\n";
	    $rnazid++;
	
	}

    }
}
