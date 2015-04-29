#!/usr/bin/perl -w
use Env qw(USER);
# tblastn -links -group
use strict;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::PrimarySeq;
use File::Spec;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::Tools::Genewise;

my %compress = ('gz'   => 'zcat',
		'bz2'  => 'bzcat',
		);
my $max_models= 100;
my $maxintron = 800;
my $informat = 'tfastx'; # wublast mformat 3 is default
my $method = 'exonerate';
my %protein2genome = ('exonerate' =>
		      'exonerate -m p2g --bestn 1 --joinfilter 1 --verbose 0 --maxintron 3000 --minintron 20 -q %s -t %s --ryo ">%%qi__%%ti (%%tab - %%tae) score=%%s rank=%%r\n%%tcs\n" --showcigar no --showvulgar no --showalignment no --refine region',
#--model protein2genome --bestn 1 --refine region  --showvulgar yes --softmaskquery yes --softmasktarget yes --minintron %d -q %s -t %s --maxintron 3000 --ryo ">%%qi length=%%ql alnlen=%%qal\n>%%ti length=%%tl alnlen=%%tal\n" --showalignment no --showtargetgff | perl /rhome/jstajich/src/genome-scripts/data_format/process_exonerate_gff3.pl --type Protein',

		      'genewise' => 'genewise -silent -quiet -para -genesf %s -u %d -v %d %s %s',
		      
		      );

my $tmpdir = "/tmp/protein2genome.$USER$$";
my ($genome,$pepfile);
my $window = 1500;  # window around target region
my $gffver = 3;
my $output;
my $debug = 0;
my $best_only = 1;
my $evalue = 1e-5;
my $once = 0;
my $basename;
GetOptions(
	   'o|out:s'         => \$output,
	   'gff|gffver:s'    => \$gffver,
	   'g|d|db|genome:s' => \$genome,
	   'p|q|qdb|pep:s'   => \$pepfile,
	   'f|format:s'      => \$informat,
	   't|tmpdir:s'      => \$tmpdir,
	   'modelmax:s'      => \$max_models,
	   'w|window:i'      => \$window,
           'b|basename:s'    => \$basename,
	   'maxintron:i'     => \$maxintron,
	   'm|method:s'      => \$method,
	   'v|debug!'        => \$debug,
           'once!'           => \$once,
	   'best!'           => \$best_only,
    );

mkdir($tmpdir) unless -d $tmpdir;
mkdir("$tmpdir/nt") unless -d "$tmpdir/nt";
mkdir("$tmpdir/pep") unless -d "$tmpdir/pep";

my $gdb = Bio::DB::Fasta->new($genome);
my $pepdb = Bio::DB::Fasta->new($pepfile);

my $out;
if($output) {
    $out = Bio::Tools::GFF->new(-gff_version => $gffver,
				-file        => ">$output");
} else {
    $out = \*STDOUT;
#    $out = Bio::Tools::GFF->new(-gff_version => $gffver);
}

for my $file ( @ARGV) {
    my $fh;
    if( $file =~ /\.(gz|bz2)$/ ) {
	open( $fh => "$compress{$1} $file |") || die $!;
    } else {
	open( $fh => $file) || die $!;
    }
    if( lc($informat) eq 'wutab') {
	my %lastseen;
	my @hsps;
	my $i = 0;
	
	while(<$fh>) {
	    chomp;
	    next if(/^\#/);
	    my ($q,$h,$e,$N,$Sprime,$S,
		$alignlen,$nident,$npos,
		$nmism, $pcident,$pcpos,$qgaps,$qgaplen,
		$sgaps,$sgaplen,
		$qframe,$qstart,$qend,
		$sframe,$sstart,$send,
		$group, $links) = split(/\t/,$_);
	    warn("group is $group links are $links for $q $h $e\n") if $debug;
	    $links =~ s/[\(\)]//g;
	    if( @hsps &&
		( $lastseen{'query'} ne $q ||
		  $lastseen{'hit'}   ne $h ) ) {
		for my $group ( keys %{$lastseen{'groups'}} ) {
		    &make_pair($lastseen{'query'},
			       $lastseen{'hit'},
			       [ map { $hsps[$_ - 1] } split(/\-/,$group)]);
		}
		@hsps = ();
		$lastseen{'groups'} = {}; # reset the groups    
		#last if $debug && $i++ > 3;
	        last if $once;
	    }

	    $lastseen{'groups'}->{$links}++;
	    ($sstart,$send) = sort { $a <=> $b} ($sstart,$send);
	    push @hsps, { qstart => $qstart,
			  qend   => $qend,
			  hstrand => substr($sframe,0,1),
			  hstart => $sstart,
			  hend   => $send};
	    $lastseen{'query'} = $q;
	    $lastseen{'hit'}   = $h;
	}
	if( @hsps ) {
	    for my $group ( keys %{$lastseen{'groups'}} ) {
		&make_pair($lastseen{'query'},
			   $lastseen{'hit'},
			   [ map { $hsps[$_] } split(/\-/,$group)]);
	    }
	}	
    } elsif ( lc($informat) eq 'tblastx') { # also for tfastx
	my %lastseen;
	my @hsps;
	my $i = 0;
	my %query_count;
	while(<$fh>) {
	    chomp;
	    next if(/^\#/);
	    my ($q,$h,$pid,$alnlen,$mismatches,$gapOpenCount,
		$qstart,$qend,
		$sstart,$send,
		$e,$bits) = split(/\t/,$_);
	    next if $e > $evalue;
	    if( @hsps &&
		( $lastseen{'query'} ne $q ||
		  $lastseen{'hit'}   ne $h ) ) {
		if( $query_count{$lastseen{'query'}}++ < $max_models) {
		 &make_pair($lastseen{'query'},
			   $lastseen{'hit'},\@hsps);
	        }
		@hsps = ();
	        last if $once;
	    }
	    my $hstrand = $sstart < $send ? 1 : -1;
	    ($sstart,$send) = sort { $a <=> $b} ($sstart,$send);
	    push @hsps, { qstart => $qstart,
			  qend   => $qend,
			  hstrand => $hstrand,
			  hstart => $sstart,
			  hend   => $send};
	    $lastseen{'query'} = $q;
	    $lastseen{'hit'}   = $h;
	}
	if( @hsps && $query_count{$lastseen{'query'}}++ < $max_models) {
	    &make_pair($lastseen{'query'},
		       $lastseen{'hit'},\@hsps);
	}
    }
}


sub make_pair {
    my ($q,$h,$hsps) = @_;
    my $qname = $q;
    my $protein = Bio::PrimarySeq->new(-seq => $pepdb->seq($q),
				       -id  => $q);

    my @srted = ( sort { $a->{hstart} <=> $b->{hstart} } @$hsps);

    my @sets;
    my $last;
    my $current_set = 0;
    for my $n ( @srted ) {
	if( $last && abs($last->{hend} - $n->{hstart}) > $maxintron) {
	    $current_set++;
	}
	push @{$sets[$current_set]}, $n;
	$last = $n;
    }
    push @{$sets[$current_set]}, $last;
    my $setct = 0;
    for my $set ( @sets ) {
	my ($min,$max) = ( $set->[0]->{hstart},
			   $set->[-1]->{hend});
	warn("set has this item ", $set->[0]->{hstart}, " ", $set->[0]->{hend}, " and ", scalar @$set, " HSPs\n") if $debug;
	my $strand = $set->[0]->{hstrand} eq '+' ? 1 : -1;
	$min -= $window; 
	$min = 1 if $min < 1;
	$max += $window;
	warn("h is $h\n") if $debug;
	my $len = $gdb->length($h);
	if( ! $len ) {
		die " cannot find seq $h in db\n";
	}
	$max = $len if $max > $len;

	my ( $chrom_frag,$run_cmd);
	
	if( $method eq 'genewise' ) {
	    $chrom_frag = $gdb->get_Seq_by_acc($h);
	} else {
	    $chrom_frag = Bio::PrimarySeq->new(-id => 
					       sprintf("%s_%d_%d_%s",
						       $h,$min,$max,
						       $set->[0]->{hstrand}),
					       -seq => $gdb->seq($h,$min,$max));
	}
	my $pepid = $protein->id;
	if( ! $pepid ) {
	 die " cannot find a pepid in seq for $q\n";
	}
	$pepid =~ s/\|/_/g;
	my $pfile = File::Spec->catfile($tmpdir,'pep',$pepid);
	warn($protein->id, " ", $chrom_frag->id, "\n") if $debug;
	if( ! -f $pfile ) {
	    Bio::SeqIO->new(-format => 'fasta',
			    -file   => ">$pfile")->write_seq($protein);
	} 
	my $fragid = $chrom_frag->id;
	$fragid =~ s/\|/_/g;
	my $gfile = File::Spec->catfile($tmpdir,'nt',$fragid);
	if( ! -f $gfile ) {
	    Bio::SeqIO->new(-format => 'fasta',
			    -file   => ">$gfile")->write_seq($chrom_frag);
	}

	if( $method eq 'genewise' ) {
	    my $strnd = $strand eq '-' ? '-trev' : '-tfor';
	    $run_cmd = sprintf($protein2genome{$method},
			       $strnd,$min,$max,
			       $pfile,$gfile);
	    
	} elsif( $method eq 'exonerate')  {
	    $run_cmd = sprintf($protein2genome{$method},
#			   $maxintron,
			       $pfile,$gfile);
	} else {
	    die "method $method unrecognized";
	}
	warn("$run_cmd ($method)\n") if $debug; 
	open( my $runfh => "$run_cmd |") || die $!;
	
	if( $method eq 'genewise' ) {
	    my $gw = Bio::Tools::Genewise->new(-fh => $runfh);
	    while (my $gene = $gw->next_prediction){
		$gene->primary_tag('gene');
		$gene->source_tag('Genewise');	
		$out->write_feature($gene);
		my @transcripts = $gene->transcripts;
		foreach my $t(@transcripts){
		    $t->primary_tag('mRNA');
		    $t->source_tag('Genewise');
		    $out->write_feature($t);
		    my @exons =  $t->exons;	    
		    foreach my $e(@exons){
			$e->remove_tag('supporting_feature');
			$e->source_tag('Genewise');
			$e->primary_tag('CDS');
			if( $e->has_tag('phase') ) {
			    $e->frame($e->get_tag_values('phase'));
			    $e->remove_tag('phase');
			}
			$out->write_feature($e);
		    }
		}
	    }
	} elsif( $method eq 'exonerate'  ) {
	    while(<$runfh>) {
		print $out $_;
	    } 
#	$out->write_feature();
	}
    }
}
