#!/usr/bin/perl -w
use strict;
use Bio::Tools::RepeatMasker;

my $file = shift;
my $in = Bio::Tools::RepeatMasker->new(-file => $file);
print "##gff-version 3\n";
my %count;
my $count = 0;
while( my $result = $in->next_result ) {
    my ($queryf,$targetf) = ($result->feature1,$result->feature2);
    my ($target,$tstart,$tend) = $queryf->get_tag_values('Target');
    print join("\t", $queryf->seq_id,
	       'repeatmasker',
	       'match_part',
	       $queryf->start,
	       $queryf->end,
	       $queryf->score,
	       $queryf->strand > 0 ? '+' : '-',
	       '.',
	       sprintf("ID=RMf.%d;Name=RM.%s.n%d;Note=%s;Target=%s+%d+%d",
		       $count++,
		       $target,$count{$target}++,
		       $queryf->primary_tag,
		       $target,$tstart,$tend)),"\n";
}
