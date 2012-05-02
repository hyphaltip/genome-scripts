#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Getopt::Long;
use Env qw(HOME);
# here is cmdline example 
# annot_xfer_mercatorblocks.pl -f gene:Broad -t gene:JGI \
#    -dt Batrachochytrium_dendrobatidis_JAM81_5 \
#    -df Batrachochytrium_dendrobatidis_JEL423_1 -m alignments \
#    -gt JAM81 -gf JEL423 > JAM81-JEL423.mercator_orthologs.tab

my ($user,$pass,$dbname_from,$dbname_to,$dbname_to2,$host);
$host ='localhost';
my $prefix;
my $debug = 0;

# move this to a YAML file?
my %genomes = ('Nc' => ['js_selker_Nc', 'gene:NC10_CALLGENES_FINAL_5'],
	       'Nt' => ['js_selker_Nt', 'gene:JGI'],
	       'Nd' => ['js_selker_Nd', 'gene:JGI']);
my @genome_order = qw(Nc Nt Nd);
my ($mercator_dir,$output);
my $ref_genome = 'Nc';
GetOptions(
	   'v|verbose!' => \$debug,
	   'u|user:s' => \$user,
	   'p|pass:s' => \$pass,
	   'host:s'   => \$host,
           'r|ref:s'  => \$ref_genome,
	   'm|mercator:s' => \$mercator_dir, # alignment dir
	   'output:s'  => \$output,
	   );

die "unknown ref genome '$ref_genome' - expected one in the ", join(",", keys %genomes), " set\n" unless $genomes{$ref_genome};

unless( defined $mercator_dir && -d $mercator_dir ) {
    die("cannot open $mercator_dir, provide with -m or --mercator\n");
}

if( $output && $output ne '-' ) { 
    open(my $fh => ">$output" ) || die $!;
    $output = $fh;
} else {
    $output = \*STDOUT;
}


($user,$pass) = &read_cnf($user,$pass) unless $pass && $user;
my %dbh;
for my $genome ( keys %genomes ) {
    my $dsn = sprintf('dbi:mysql:database=%s;host=%s',$genomes{$genome}->[0],$host);
    $dbh{$genome} = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::mysql',
						    -dsn     => $dsn,
						    -user    => $user,
						    -password => $pass,
	);
    
}
if( ! defined $dbh{$ref_genome} ) {
    die("unknow refgenome $ref_genome\n");
}
my $ref_dbh = $dbh{$ref_genome};
print $output join("\t", ( map { my $sp = $_;
				 map { sprintf("%s_%s",$_, $sp) } 
				 qw(GENE MRNA CHROM START STOP STRAND)  }
			   @genome_order),'SINGLE_MATCH'), "\n";

my $iter = $ref_dbh->get_seq_stream(-type => $genomes{$ref_genome}->[1]);

while( my $gene = $iter->next_seq ) {

    my $name = $gene->name;
    my ($mRNA) = $gene->get_SeqFeatures('mRNA'); # 1st mRNA for now
    my $t_name;
    if( ! defined $mRNA ) {
	warn("no mRNA for $name\n");
	$t_name = $name;
    } else {
	$t_name = $mRNA->name;
    }
    # This program requires the sliceAlignment program from
    # MERCATOR which is part of Colin Dewey's tools
    my $arg = sprintf("sliceAlignment %s %s %s %d %d %s",
		      $mercator_dir, $ref_genome,
		      $gene->seq_id, $gene->start, $gene->end,
		      $gene->strand > 0 ? '+' : '-');
    warn("arg is $arg\n") if $debug;
    open(my $map => "$arg 2>/dev/null | grep '>' | " ) || die "Cannot open slice with $arg\n";
    my $seen = 0;
    my %segments;
    while(<$map>) {
	# parsing output from sliceAlignment, all we really want is the 
        # matching interval in the other genome (genome_to is the prefix)
	if( /^>(\S+)\s+([^:]+):(\d+)\-(\d+)([+-])/ ) {
	    my ($genome,$chrom,$start,$end,$strand) = ($1,$2,$3,$4,$5);	
	    if( $genome ne $ref_genome ) {
		my $dbh = $dbh{$genome};
		$seen = 1;
		# get that segment in the 'TO' genome
		my ($segment) = $dbh->segment($chrom,$start,$end);
		# extract the gene(s) in this interval
		if( ! defined $segment ) { 
			warn("cannot find segment in the DB for $chrom:$start..$end ($genome)\n");
			next;
		}
		my @genes = $segment->features(-type => $genomes{$genome}->[1]);
		if( @genes ) {
		    for my $g ( @genes ) {
			my ($to_mRNA) = $g->get_SeqFeatures('mRNA');
			# print out the genes that fall in syntenic interval
			push @{$segments{$genome}}, [$g->name,$to_mRNA->name,
						     $g->seq_id,$g->start,$g->end,
						     $g->strand];
		    }
		} else {
		    push @{$segments{$genome}}, ['','','','','','NO_GENES_IN_INTERVAL'];
		}
	    }
	}
    }
    
    # merge positions
    my %merged;
    my $single = 'yes';
    for my $gn ( keys %segments ) {
	my (@gnames,@mnames,@strands,$multi);
	my ($seqid,$min,$max);
	for my $dat ( @{$segments{$gn}} ) {
	    if( $dat->[5] !~ /^NO_/ ) {
		push @gnames, $dat->[0];
		push @mnames, $dat->[1];
		if( defined $seqid && $dat->[2] ne $seqid) {
		    warn("mixed chromosomes for an interval??? @gnames | @mnames\n");
		}
		$seqid = $dat->[2];		
		$min = $dat->[3] if (! defined $min || $dat->[3] < $min);
		$max = $dat->[4] if (! defined $max || $dat->[4] > $max);
	    }
	    push @strands, $dat->[5];		
	}
	$merged{$gn} = [ join(",", @gnames),
			 join(",", @mnames),
			 $seqid || '',  $min || '', $max || '', 
			 join(",", @strands)];
	$single = 'no' if @gnames > 1;
    }
    print $output join("\t",
		       $name, $t_name,$gene->seq_id,
		       $gene->start,$gene->end,
		       $gene->strand,
		       ( map { @{$merged{$_} || ['','','','','','NO_SYNTENIC_ALIGNMENT']} } grep { $_ ne $ref_genome } @genome_order),
		       $single),"\n";
    last if $debug;
}


# read the .my.cnf file
sub read_cnf {
    my ($user,$pass) = @_;
    if( -f "$HOME/.my.cnf") {
        open(IN,"$HOME/.my.cnf");
        while(<IN>) {
            if(/user(name)?\s*=\s*(\S+)/ ) {
                $user = $2;
            } elsif( /pass(word)\s*=\s*(\S+)/ ) {
                $pass = $2;
            }
        }
        close(IN);
    }
    return ($user,$pass);
}
