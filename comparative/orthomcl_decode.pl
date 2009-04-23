#!/usr/bin/perl -w
# Allows re-running of MCL from a starting input .mtx file to produce a new
# mcl with different inflation parameters etc
# 
# ie
# cd MYorthomcl_output
# cd tmp
# mcl -I 2.2 -te 6 all_ortho.mtx -o all_ortho.I_22.mcl
# orthomcl_decode.pl file.gg all_ortho.I_22.mcl all_ortho.idx > all_orthomcl.I_22.out

use strict;
use Getopt::Long;

my $gg_file     = shift || die;
my $mcl_file    = shift || die;
my $mcl_index   = shift || die;
my (%gindex,%gindex2,%mcl_backindex,@mcl_index);

# read usr_gg_file
open (GG,$gg_file) or dieWithUnexpectedError("can't open $gg_file!");
while (<GG>) {
    s/[\r\n]//g;
    if (/(\S+)\(\d+\):(.*)/) {
	my ($taxon) = ($1);
	my @genes=split (" ",$2);
	push (@{$gindex{$taxon}},@genes);
	foreach (@genes) {$gindex2{$_}=$taxon;}
	#$totalgeneno+=scalar(@{$gindex{$taxon}});
	warn($_);
    } elsif (/(\S+):(.*)/) {
	my $taxon=$1;
	my @genes=split (" ",$2);
	push (@{$gindex{$taxon}},@genes);
	for (@genes) {$gindex2{$_} = $taxon;}
	#$totalgeneno+=scalar(@{$gindex{$taxon}});
    }
}

# PARSE MCL INDX FILE
{ 
    open(MCLIDX, $mcl_index) || die $!;
    while(<MCLIDX>) {
	my ($id,$gene)        = split;
	$mcl_index[$id]       = $gene;
	$mcl_backindex{$gene} = $id;
    }    
}

# WRITE ORTHOMCL out
{
    open (MCL,$mcl_file) or dieWithUnexpectedError("can't open $mcl_file");
    my $last=0;
    my $lastline='';
    my @mcl;
    while (<MCL>) {
#		chomp;chop;  #bug, reported by Robson Francisco de Souza.
	#this may result in wrong gene index by chopping the last digit of line
	#for new format of MCL output. Not a bug for older versions, e.g. mcl-02-063
	#replaced to the following line: substitute the '$' with nothing, and 
	#end reading data in the end of data block, preventing loading comment
	#rows.
	chomp; s/\$$//; last if (/^\)/ && scalar(@mcl));	
	if (/^(\d+)\s+(.*)\$/) {
	    $mcl[$last]=$lastline;
	    $last=$1;$lastline=$2;
	}
	elsif (/^(\d+)\s+(.*)/) {
	    $mcl[$last]=$lastline;
	    $last=$1;$lastline=$2;
	}
	elsif (/^\s+/) {$lastline.=$_;}
    }
    $mcl[$last]=$lastline;
    close (MCL);

    my $orthomcl_cluster_id = 1;
    foreach my $mcl_cluster_id (0..$last) {
	$mcl[$mcl_cluster_id]=~s/\s+/ /g;
	my @g=split (" ",$mcl[$mcl_cluster_id]);
	next unless (scalar(@g)>=2);
	my @taxa=();
	foreach (@g) {
	    my $taxon = $gindex2{$mcl_index[$_]};
	    my $presence=0;
	    foreach (@taxa) {
		if ($_ eq $taxon) {$presence=1;}
	    }
	    push (@taxa,$taxon) unless ($presence);
	}

	print "ORTHOMCL".$mcl_cluster_id."(".scalar(@g)." genes,".scalar(@taxa)." taxa):	";
	foreach (@g) {
	    print " $mcl_index[$_]($gindex2{$mcl_index[$_]})";
	}
	print "\n";
#No Species Cutoff
    }

}

