#!/usr/bin/perl -w
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

my $maxintron = 500;
my $informat = 'wutab'; # wublast mformat 3 is default
my $method = 'genewise';
my %protein2genome = ('exonerate' =>
		      "exonerate --model protein2genome --bestn 1 -E yes --showtargetgff  --maxintron %d -q %s -t %s ",
		      'genewise' => 'genewise -silent -quiet -para -genesf %s -u %d -v %d %s %s',
		      
		      );

my $tmpdir = '/tmp/protein2genome';
my ($genome,$pepfile);
my $window = 500;  # window around target region
my $gffver = 3;
my $output;
GetOptions(
	   'o|out:s'         => \$output,
	   'gff|gffver:s'    => \$gffver,
	   'g|d|db|genome:s' => \$genome,
	   'p|q|qdb|pep:s'   => \$pepfile,
	   'f|format:s'      => \$informat,
	   't|tmpdir:s'      => \$tmpdir,
	   'w|window:i'      => \$window,
	   'maxintron:i'     => \$maxintron,
	   'm|method:s'      => \$method,
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
    $out = Bio::Tools::GFF->new(-gff_version => $gffver);
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
	    next if(/^#/);
	    my ($q,$h,$e,$N,$Sprime,$S,
		$alignlen,$nident,$npos,
		$nmism, $pcident,$pcpos,$qgaps,$qgaplen,
		$sgaps,$sgaplen,
		$qframe,$qstart,$qend,
		$sframe,$sstart,$send,
		$group,$links) = split(/\t/,$_);

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
#		last if $i++ > 3;
	    }
	    $links =~ s/[\(\)]//g;
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
    }
}


sub make_pair {
    my ($q,$h,$hsps) = @_;
    my $qname = $q;
    my $protein = Bio::PrimarySeq->new(-seq => $pepdb->seq($q),
				       -id  => $q);
    my @srted = ( sort { $a->{hstart} <=> $b->{hstart} } @$hsps);
    my ($min,$max) = ( $srted[0]->{hstart},
		       $srted[-1]->{hend});
    my $strand = $srted[0]->{hstrand} eq '+' ? 1 : -1;
    $min -= $window; 
    $min = 1 if $min < 1;
    $max += $window;
    my $len = $gdb->length($h);
    $max = $len if $max > $len;

    my ( $chrom_frag,$run_cmd);
    
    if( $method eq 'genewise' ) {
	$chrom_frag = $gdb->get_Seq_by_acc($h);
    } else {
	$chrom_frag = Bio::PrimarySeq->new(-id => 
			       sprintf("%s_%d_%d_%s",
				       $h,$min,$max,
				       $srted[0]->{hstrand}),
			       -seq => $gdb->seq($h,$min,$max));
    }
    my $pfile = File::Spec->catfile($tmpdir,'pep',$protein->id);
    #warn($protein->id, " ", $chrom_frag->id, "\n");
    if( ! -f $pfile ) {
	Bio::SeqIO->new(-format => 'fasta',
			-file   => ">$pfile")->write_seq($protein);
    }
    my $gfile = File::Spec->catfile($tmpdir,'nt',$chrom_frag->id);
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
			   $maxintron,
			   $pfile,$gfile);
    } else {
	die "method $method unrecognized";
    }

    open( my $runfh => "$run_cmd |") || die $!;
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
		$out->write_feature($e);
	    }
	}
    }

}
