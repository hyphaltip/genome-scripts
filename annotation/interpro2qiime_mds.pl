#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use GO::Parser;
use GO::AppHandle;
use Bio::DB::Taxonomy;

my %altids = ('GO:0006467' => 'GO:0003756');

my ($host,$user,$pass) = ('puffball.fungalgenomes.org',
			  'GOuser',
			  'G0user');
my $taxdb = Bio::DB::Taxonomy->new(-source => 'entrez');    

#my $parser = GO::Parser->new({handler => 'obj',use_cache=>1});

my $iprin;
my $outdir= ".";
my $ext = 'IPROUT.tsv';
my $gofile;
GetOptions(
    'o|out:s' => \$outdir,
    'i|ipr:s' => \$iprin,
    'ext:s'   => \$ext,
    'go:s'    => \$gofile,
    );

die "need an iprfile " unless( defined $iprin);

my $apph = GO::AppHandle->connect(-dbname => 'GO',
				  -dbhost => $host,
				  -dbuser => $user,
				  -dbauth => $pass);

my %cache;
my (%gene2pathways,%gene2domain);
my $stem;
unless ($iprin =~ /(\S+)\.$ext$/) {
    warn("extension $ext does not exist for file $iprin\n");
    $stem = $iprin;
} else {
    $stem = $1;
}
open(my $fh => $iprin ) || die "cannot open file: $iprin!";
while(<$fh>) {
    my ($seqid, $md5, $length, $db, $dbacc, $dbdesc, $alnstart,$alnend,
	$evalue, $status, $date,$iprdomain, $iprdesc,
	$go, $pathways) = split(/\t/,$_);
    $gene2domain{$db}->{$seqid}->{$dbacc} = $dbdesc;
    
    if( $go ) {
	for my $goid ( split(/\|/,$go) ) {
	    $gene2domain{GO}->{$seqid}->{$goid}++;
	}
    }
    if( $pathways ){
	for my $pathway ( split(/\|/,$pathways) ) {
	    if( $pathway =~ /(\S+):\s+(\S+)/ ) {
		$gene2pathways{$seqid}->{$1}->{$2}++;
	    }
	}
    }
}
open(my $ofh => ">$outdir/$stem\_go_pathway.txt") || die $!;
for my $db ( qw(GO) ) {
    for my $seqid ( sort keys %{$gene2domain{$db}} ) {
	my @row;
	for my $go ( sort keys %{$gene2domain{$db}->{$seqid}} ) {	    
#	    my $graph = $apph->get_node_graph($go, 4);
#		my $iter = $graph->create_iterator;
	    if( ! exists $cache{$go} ) {
		$go = $altids{$go} if( exists $altids{$go} );
		my $term = $apph->get_term({acc=>$go});
		if( ! $term ) {
		    warn("cannot find term $go\n");
		    next;
		}
		if( $term->get_code_from_namespace eq 'C' ) {
			next;
		}
		my @paths;
		for my $path (@{$apph->get_paths_to_top({acc=>$go})}) {
		    my @all = map { $_->name } reverse @{$path->term_list};
		    my @path = splice(@all,1,5);
			# need to order these later
#		while (my $ni = $iter->next_node_instance) {
#		    my $depth = $ni->depth;
#		    my $term = $ni->term;
		    #my $reltype = $ni->parent_rel->type;
#		    push @path, $term->name;	    		
		    push @paths, [@path];
		}
		my ($longest) = sort { scalar @{$b} <=> scalar @{$a}} @paths;
		$cache{$go} = $longest;
	    }
	    push @row, join("; ", @{$cache{$go}},$go);
	}
	unless( @row ) {
	    push @row, "Unknown";
	}
	print $ofh join("\t", $seqid, @row),"\n";
    }
}
