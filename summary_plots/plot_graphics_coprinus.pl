#!/usr/bin/perl -w
# $Id$

use strict;
use File::Spec;
use Bio::SeqIO;
use List::Util qw(sum max);
use Cwd;
use Getopt::Long;
use Bio::Root::RootI;
use GD qw(gdGiantFont gdLargeFont gdMediumBoldFont gdSmallFont gdTinyFont);
use SVG;

# SHOULD WE GENERATE THE PLOT WITH GD OR SVG?
# This should just be a command line parameter

my $GD  = 0;

mkdir("png") unless -d "png";
mkdir("svg") unless -d "svg";

use constant TICK_SCALE => 1_000_000; # draw tickmark every 1M bp

my $DIR = 'plot';
my $dbfile;

GetOptions(
	   'd|dbfile:s' => \$dbfile,
	   'dir:s'      => \$DIR,
	   'gd!'        => \$GD,
	   );

if( ! $dbfile ) {
    $dbfile = File::Spec->catfile($DIR, 'genome.fa');
}


# recalculate values from all files...
my $USE_CACHED = 0;


# DATA SOURCES
use constant GENES    => 'coprinus_gene_summary.tab';
use constant ORTHOS   => 'coprinus_orthologs.tab';
use constant ORPHANS  => 'coprinus_orphans.tab';
use constant PARALOGS => 'coprinus_paralogs.tab';
use constant SSGENES  => 'coprinus_only.tab';

use constant FISH_SIGNIF_BLOCKS => 'coprinus_FISH_synteny.tab';
use constant FISH_BLOCKS_ALL    => 'coprinus_FISH_synteny_original.tab';
use constant LB_BLOCKS     => 'cc-lb_synteny.tab';
use constant PC_BLOCKS     => 'cc-pc_synteny.tab';

use constant REPEATS  => 'coprinus_rptmask.tab';
use constant INTERGENIC => 'coprinus_intergenic_summary.tab';

use constant RECOMBINATION_RATES => 'recombination_rates.tab';
use constant SSR => 'ssr.tab';
use constant TRNA => 'coprinus_trna.tab';
use constant TE_REPEAT => 'coprinus_TE_repeats.tab';

use constant CEN_TEL => 'coprinus_CEN_TEL.tab';

use constant DS_DN_AVG => 'coprinus_avg_ds.tab';
use constant DS_DN_AVG_INVERT => 'coprinus_avg_ds_invert.tab';

# TRACK COLORS - COLLECT HERE FOR EASY ACCESS
use constant BLOCKS_COLOR  => 'steelblue';
use constant BLOCKS_COLOR_REV  => 'goldenrod';

use constant RECOMBINATION_HOT_COLOR  => 'red';
use constant RECOMBINATION_COLD_COLOR  => 'steelblue';
use constant RECOMBINATION_NEUTRAL_COLOR  => 'darkgrey';

use constant GENES_COLOR     => 'darkorange';
use constant ORPHANS_COLOR   => 'orange';
use constant ORTHOS_COLOR    => 'varP2';
use constant PARALOGS_COLOR  => 'red';

use constant KINASE_DOMAINS_COLOR   => 'steelblue';
use constant P450_DOMAINS_COLOR     => 'steelblue';
use constant WD40_DOMAINS_COLOR     => 'steelblue';

use constant INTERGENIC_COLOR => 'varP4';
use constant TRNA_COLOR      => 'limegreen';
use constant TE_REPEAT_COLOR => 'brown';
use constant REPEATS_COLOR   => 'varB1';

use constant FISH_BLOCKS_COLOR => 'lightgreen';
use constant FISH_SIGNIFBLOCKS_COLOR => 'darkgreen';

use constant CENTROMERE_COLOR => 'black';
use constant TELOMERE_COLOR   => 'red';

use constant DS_COLOR => 'blue';
use constant DN_COLOR => 'red';
# Font sizes...
use constant LABEL_SIZE    => '12';

# These are not being used right now
my %labels = ( 
	       genes           => 'Genes / 50 kb',
	       intergenic      => 'Intergenic distances / 50 kb',
	       gmap            => 'gmap vs pmap',
	       orthologs       => 'Orthologous genes / genes / 50 kb',
	       orphans         => 'Orphan genes / genes / 50 kb',
	       paralogs        => 'Paralogous genes / genes / 50 kb',
	       mercator_blocks => 'Syntenic blocks',
	       ortholog_blocks => 'Ortholog blocks',
	       pc_blocks       => 'Pchr Mercator Syntenic blocks',
	       lb_blocks       => 'Lbic Mercator Syntenic blocks',
	       repeats         => 'All Repetitive elements / 50 kb',
	       te_repeat       => 'Transposable elements / 50 kb',
	       recombination_rates => 'Recombination rate',
	       ssr             => 'SSRs',
	       dS              => 'Paralogs dS',
	       dSinv              => 'Paralogs 1/dS',
	       dN              => 'Paralogs dN',
	       tRNA            => 'tRNA genes / 50 kb',
	       fish_blocks     => 'FISH synteny blocks',
	       fish_blocks_all  => 'FISH synteny blocks',
	       centromere      => 'Telomere and Centromere location',
	       );
my %ylabels = (
	       genes           => 'genes',
	       intergenic      => 'intergenic',
	       orthologs       => 'ortholog-genes',
	       orphans         => 'orphan-genes',
	       paralogs        => 'paralog-genes',
	       ssgenes         => 'species-specific-genes',
	       ortholog_blocks => 'Ortholog blocks',
	       pc_blocks       => 'Pchr Syntenic blocks',
	       lb_blocks       => 'Lbic Syntenic blocks',
	       repeats         => 'repetitive elements-50kb',
	       te_repeat       => 'te-repeats',
	       kinase_domains  => 'kinasedomains',
	       p450_domains    => 'p450domains',
	       wd40_domains    => 'wd40domains',
	       recombination_rates       => 'recombination',
	       ssr             => 'SSRs',
	       dS              => 'dS',
	       dSinv              => 'dS-inv',
	       dN              => 'dN',
	       tRNA            => 'tRNA',
	       fish_blocks     => 'fish_blocks',
	       fish_blocks_all => 'fish_blocks_all',
	       centromere => 'centromere',
	       );

my %units = (
	     genes           => 'genes/50kb',
	     intergenic      => 'Intergenic distances/50kb',
	     gmap            => 'gmap vs pmap',
	     orthologs       => 'orthologs/genes/50kb',
	     orphans         => 'orphans/genes/50kb',
	     paralogs        => 'paralogs/genes/50kb',
	     ssgenes         => 'species-specific/genes/50kb',
	     ortholog_blocks => 'Ortholog blocks',
	     mercator_blocks => 'Syntenic blocks',
	     pc_blocks       => 'Pchr Mercator Syntenic blocks',
	     lb_blocks       => 'Lbic Mercator Syntenic blocks',
	     repeats         => 'All repetitive elements/50kb',
	     te_repeat      => 'TE repeats / 50kb',
	     kinase_domains  => 'Kinase Domains',
	     p450_domains    => 'P450 Domains',
	     wd40_domains    => 'WD40 Domains',
	     recombination_rates       => 'Recombination rates',
	     ssr             => 'SSRs',
	     dS              => 'dS/50kb',
	     dSinv              => '1/dS/50kb',
	     dN              => 'dN/50kb',
	     tRNA            => 'tRNA genes',
	     fish_blocks     => 'FISH blocks',
	     fish_blocks_all     => 'FISH blocks',
	     centromere      => 'Telomere and Centromere location',
	     );

# IMAGE CONSTANTS
use constant TOP           => 40;
use constant TRACK_HEIGHT  => 40;
use constant TRACK_SPACE   => 40;
use constant TOTAL_TRACKS  => 11.5;

# OLD VALUES FOR TOP ALIGNED LABEL
use constant TRACK_LEFT    => 35;
use constant WIDTH         => 800;

# Values for left aligned labels
#use constant WIDTH         => 1300;
#use constant TRACK_LEFT    => 200;

# Absolute chromosome lengths

my %CHROMS;
{
    my $in = Bio::SeqIO->new(-format => 'fasta',-file=>$dbfile);
    my $i = 0;
    while( my $seq = $in->next_seq ) {
	my $id = $seq->display_id;
	$CHROMS{$id}= [$i++, $seq->length];
    }
}

my ($largest_chrom) = sort { $b<=>$a } map {$_->[1]} values %CHROMS;  # inefficient sort

my $CONFIG = {};

# Create a bunch of objects containing the parsed data...
my $genes = Windows->new('genes');
$genes->parse( File::Spec->catfile($DIR,GENES),
		  qw(chrom chrom_start));

my $repeats = Windows->new('repeats');
$repeats->parse(File::Spec->catfile($DIR,REPEATS),
		qw(chrom start));

my $te = Windows->new('te_repeat');
$te->parse(File::Spec->catfile($DIR,TE_REPEAT),
	   qw(chrom start));

my $orthos = Windows->new('orthologs');
$orthos->parse(File::Spec->catfile($DIR,ORTHOS),
	       qw(src src_start));
$orthos->normalize($genes);

my $orphans = Windows->new('orphans'); # orphans
$orphans->parse(File::Spec->catfile($DIR,ORPHANS),
		qw(chrom chrom_start));
$orphans->normalize($genes);

my $paralogs = Windows->new('paralogs'); #paralogs
$paralogs->parse(File::Spec->catfile($DIR,PARALOGS),
		qw(chrom chrom_start));
$paralogs->normalize($genes);

# my $ssgenes = Windows->new('ssgenes'); #species-specificgenes
# $ssgenes->parse(File::Spec->catfile($DIR,SSGENES),
# 		qw(chrom chrom_start));
# $ssgenes->normalize($genes);

my $intergenic = Windows->new('intergenic'); # intergenic distances
$intergenic->parse(File::Spec->catfile($DIR,INTERGENIC),
		   qw(chrom chrom_start length), sub { $_[0] < 1000 && $_[0] > 10 } );

my $recombination_rates = Windows->new('recombination_rates'); # recombination rates
$recombination_rates->parse_recombination_rates(File::Spec->catfile
						($DIR,RECOMBINATION_RATES));

my $ssrs = Windows->new('ssr'); #SSRs
$ssrs->parse_range(File::Spec->catfile($DIR,SSR),
		   qw(chrom start stop));

my $fish_blocks_signif = Windows->new('fish_blocks'); #FISH blocks
$fish_blocks_signif->parse_blocks(File::Spec->catfile($DIR,FISH_SIGNIF_BLOCKS),
			   qw(chrom start stop block));
my $fish_blocks_all = Windows->new('fish_blocks_all'); #FISH blocks
$fish_blocks_all->parse_blocks(File::Spec->catfile($DIR,FISH_BLOCKS_ALL),
			   qw(chrom start stop block));

my $trna = Windows->new('tRNA'); # tRNAs
$trna->parse(File::Spec->catfile($DIR,TRNA),
	     qw(chrom start));
# Just make the settings has global so I don't have to worry about
# passing it around 'cuz that just sucks

my $centromere = Windows->new('centromere'); # CENTROMERE & TELOMERE locations
$centromere->parse_blocks(File::Spec->catfile($DIR,CEN_TEL),
			  qw(chrom start stop type));

#my $dS = Windows->new('dS'); # dS
#$dS->parse(File::Spec->catfile($DIR,DS_DN_AVG),
#	   qw(scaffold start_position dS));

my $dS_invert = Windows->new('dSinv'); # dS
$dS_invert->parse(File::Spec->catfile($DIR,DS_DN_AVG_INVERT),
	   qw(scaffold start_position dS_inv));

#my $dN = Windows->new('dN'); # dS
#$dN->parse(File::Spec->catfile($DIR,DS_DN_AVG),
#	   qw(scaffold start_position dN));

my $settings;

for my $chrom (sort { $CHROMS{$a}->[0] <=> $CHROMS{$b}->[0] } 
	       keys %CHROMS) {
    next if $chrom =~ /^U/;
    warn( "Generating $chrom plot...\n");
    # The class of the img object depends on what we are trying to print.
    my $img = establish_image();
    $settings = establish_settings($img);
    # Make all images the same size - scaled to the largest chromosome.
    # Calculate the yscale
    my $xscale = (WIDTH - (TRACK_LEFT * 2))/$largest_chrom;

    $CONFIG = {};
    $CONFIG = { chrom  => $chrom,
		xscale => $xscale,
		count  => 0,
		img    => $img };

#    plot('repeats',$repeats,'total',REPEATS_COLOR);
    plot_centromere_telomere($centromere,CENTROMERE_COLOR, TELOMERE_COLOR);

    plot('te_repeat', $te,'total',TE_REPEAT_COLOR);
    plot('tRNA',$trna,'total',TRNA_COLOR);
    plot_recombination_rates($recombination_rates,
			  RECOMBINATION_HOT_COLOR,RECOMBINATION_COLD_COLOR,
			  RECOMBINATION_NEUTRAL_COLOR,
			  $ssrs);  
    plot('genes',$genes,'total',GENES_COLOR);
    plot('orphans',$orphans,'normalized',ORPHANS_COLOR);
    plot('orthologs',$orthos,'normalized',ORTHOS_COLOR);
    plot('paralogs',$paralogs,'normalized',PARALOGS_COLOR);
    #plot_numeric('dN',$dN,1,DN_COLOR);
    #plot_numeric('dS',$dS,3,DS_COLOR);
    plot_numeric('dSinv',$dS_invert,1,DS_COLOR);

    plot_fish_blocks($fish_blocks_all,$fish_blocks_signif,
		     FISH_BLOCKS_COLOR,
		     FISH_SIGNIFBLOCKS_COLOR);
#    plot_barplot('intergenic',$intergenic,'total',INTERGENIC_COLOR);

    # Draw some header information and the xscale, which is
    # always in megabases
    my $width   = $CHROMS{$chrom}->[1];
    my $scaled = TRACK_LEFT + ($xscale * $width);
    my $height = TOTAL_TRACKS * (TRACK_HEIGHT + TRACK_SPACE) + 15;
    my $header = "Chromosome $chrom";
    my $footer = "Megabase pairs";

    if ($GD) {
	# Place the header on the right side of the image
	$img->string(gdLargeFont,
		     $scaled - (length($header) * gdLargeFont->width) - 2,
		     5,$header,$settings->{black});
	# ...or place it on the left edge...
	#$img->string(gdGiantFont,TRACK_LEFT,
	#	       5,$header,$settings->{black});

#	$img->string(gdMediumBoldFont,$scaled/2-(length($footer)/2),
#		     $height - (gdMediumBoldFont->height) - 2,$footer,$settings->{black});
	open OUT,">png/$chrom.png";
	print OUT $img->png;
    } else {
	$img->text(
		   style=> {
		       'font' => 'Arial',
		       'font-size'  => 40,
		       'font-style' => 'bold'
		       },
		   id=>"$header",
		   x=>$scaled - (length($header) * 10) - 2,
		   y=>40)->cdata($header);

#	$img->text(
#		   style=> {
#		       'font' => 'Arial',
#		       'font-size' => 14,
#		   },
#		   id=>"Megabase pairs",
#		   x=>$scaled/2-(length($footer)/2),
#		   y=>$height)->cdata($footer);
	open OUT,">svg/$chrom.svg";
	my $out = $img->xmlify;
	print OUT $out;

    }
    close OUT;
}

sub plot {
  my ($label,$obj,$flag,$icolor,$ymax) = @_;
  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count};

  my $color = $settings->{$icolor};
  my $black = $settings->{black};
  
  my ($yscale);
  if ($ymax) {
      $yscale = TRACK_HEIGHT / $ymax;
  } else {
      $yscale = $obj->calc_y_scale($flag,TRACK_HEIGHT);
  }

  my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;

  # Create a grouping for this element for SVG-generated images

  my $bins = $obj->fetch_bins($chrom);
  for my $bin (@$bins) {
    my $total = $obj->$flag($chrom,$bin);
    my $left  = ($bin * $xscale) + TRACK_LEFT;
    my $top   = $track_baseline + (TRACK_HEIGHT - ($total * $yscale));
    my $width = 2;

    $top = ($top < $track_baseline) ? $track_baseline : $top;
    my $bottom = $track_baseline + TRACK_HEIGHT;
    if ($GD) {
      $img->rectangle($left,$top,$left+$width,$bottom,$color);
  } else {
	$img->rectangle(x=>$left,y=>$top,
			width  => $width,
			height => $bottom - $top,
			id     =>"$label-$bin",
			stroke => $color,
			fill   => $color);
    }

    # print STDERR join("\t",$left,$top,$right,$bottom),"\n";
  }
  # Draw rectangles for each
  draw_yticks($track_baseline,$yscale,$yscale,'5',$obj->{ylabel},2);
  draw_bounding($track_baseline,$obj);

  $CONFIG->{count}++;
  return;
}

# Draw y-scale ticks...
# This should just take a ymax, then calculate the ticks (evenly rounded) from that.
sub draw_yticks {
  my ($track_baseline,$yscale_top,$yscale_bot,
      $total_ticks,
      $ylabel,$sigfigs_count,$pos_neg) = @_;
  
  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count};
  my $black = $settings->{black};
  $sigfigs_count ||= 1;
  
  # Draw y-scale ticks...
  my $tick_left  = TRACK_LEFT;
  my $tick_right = TRACK_LEFT - 4;
  # HACK HACK BEGIN
  # Coprinus specific hack
  $total_ticks = 1;
  $sigfigs_count = 0;
  # HACK HACK END

# It's printing the values in backwards order
  my $interval = TRACK_HEIGHT / $total_ticks;
  my $sigfigs = "%.".$sigfigs_count."f";  # how many significant figs for the ylabel
                    # to avoid the ylabel from not being specific
  if( $yscale_top == 0 || ! defined $yscale_top ) {
      warn("yscale top is '$yscale_top'\n");
      warn(Bio::Root::RootI->stack_trace_dump());
      die;
  }
  if ( 0.9 >= ( $interval * $total_ticks) / $yscale_top ) {
      $sigfigs = "%.2f";
  } 
  if( $sigfigs_count == 0 ) {
      $sigfigs = "%2d";
  }
  
  for (my $i=0; $i <= $total_ticks;$i++) {
    my $top = $track_baseline + TRACK_HEIGHT - ($i * $interval);
    my $label = ($i * $interval) / $yscale_top;
    # do some conditional aspect of formatting
#    my $formatted_label = sprintf(".1f",$label);
    my $formatted_label = sprintf($sigfigs,$label);
    $label ||= '0';
    if ($GD) {
      $img->line($tick_left,$top,$tick_right,$top,$settings->{black});
      $img->string(gdSmallFont,
		   $tick_left- ( length($formatted_label) * gdSmallFont->width) - 4,
		   ($top - gdSmallFont->height),$formatted_label,$settings->{black});
    } else {
      $img->line(x1 => $tick_left, y1 => $top,
		 x2 => $tick_right,y2 => $top,
		 id => "$ylabel-$formatted_label ytick".$i,
		 stroke => $settings->{black},
		 fill   => $settings->{black});

      $img->text(
		 style=> {
			  'font' => 'Arial',
			  'font-size' => LABEL_SIZE,
		      },
		 id=>"$ylabel-$formatted_label-".$i,
		 x=>$tick_left-length($formatted_label) - 20, # -22
		 y=>$top+3)->cdata($formatted_label);
    }
  }
  if( $pos_neg ) {
      for (my $i=1; $i <= $total_ticks;$i++) {
	  my $top = $track_baseline + TRACK_HEIGHT + ($i * $interval);
	  my $label = ($i * $interval) / $yscale_bot;
	  # do some conditional aspect of formatting
	  my $formatted_label = sprintf("-".$sigfigs,$label);
	  $label ||= '0';

	  if ($GD) {
	      $img->line($tick_left,$top,$tick_right,$top,$settings->{black});
	      $img->string(gdSmallFont,
			   $tick_left- ( length($formatted_label) * gdSmallFont->width) - 4,
			   ($top - gdSmallFont->height),$formatted_label,$settings->{black});
	  } else {
	      $img->line(x1 => $tick_left, y1 => $top,
			 x2 => $tick_right,y2 => $top,
			 id => "$ylabel-$formatted_label ytick",
			 stroke => $settings->{black},
			 fill   => $settings->{black});

	      $img->text(
			 style=> {
			     'font' => 'Arial',
			     'font-size' => LABEL_SIZE,
			 },
			 id=>"$ylabel-$formatted_label",
			 x=>$tick_left-length($formatted_label) - 22,
			 y=>$top+3)->cdata($formatted_label);
	  }
      }
  }
}

sub plot_barplot {
    my ($label,$obj,$flag,$icolor,$ymax) = @_;
    $flag = 'value'; # hardcoded for now
    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};

    my $PLOT_BY = 'total';
    my $color = $settings->{$icolor};
    my $black = $settings->{black};


    my ($yscale);
    if ($ymax) {
	$yscale = TRACK_HEIGHT / $ymax;
    } else {
	$yscale = $obj->calc_y_scale($flag."_".$chrom,TRACK_HEIGHT);
    }

    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE)
				* $count) + TRACK_SPACE / 2;

    for my $point (@{$obj->{nonsliding}->{$chrom}}) {
	my ($pos,$val) = @{$point};

	my $left  = ($pos * $xscale) + TRACK_LEFT;
	my $top   = ($track_baseline + TRACK_HEIGHT) - ($val * $yscale);
	if ( $top < $track_baseline ) {
	    $top = $track_baseline;
	}
	if ($GD) {
	    $img->line($left,$top,$left,
		       $track_baseline + TRACK_HEIGHT,$color);
	    } else {
		$img->line(x1=>$left, y1=>$top,
			   x2=>$left, y2=>$track_baseline + TRACK_HEIGHT,
			   id=>"$label-$pos-fwd",
			   stroke=>$color,
			   fill=>$color);
	    }
	# print STDERR join("\t",$left,$top),"\n";
    }


    draw_yticks($track_baseline,
		$yscale,
		$yscale,
		'5',$obj->{ylabel},
		0, # sig figures
		0, # two panel
		);

    draw_bounding($track_baseline,$obj);
    $CONFIG->{count}++;
    return;
}

sub plot_double_barplot {
    my ($label,$obj_top,$obj_bot,
	$flag,$ifcolor,$ircolor,$ymax_top, $ymax_bot) = @_;
    $flag = 'value'; # hardcoded for now
    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};

    my $PLOT_BY = 'total';
    my $fcolor = $settings->{$ifcolor};
    my $rcolor = $settings->{$ircolor};
    my $black = $settings->{black};


    my ($yscale_top,$yscale_bot);
    if ($ymax_top) {
	$yscale_top = TRACK_HEIGHT / $ymax_top;
    } else {
	$yscale_top = $obj_top->calc_y_scale($flag."_".$chrom,TRACK_HEIGHT);
    }
    if( $ymax_bot ) {
	$yscale_bot = TRACK_HEIGHT / $ymax_bot;
    } else {
	$yscale_bot = $obj_bot->calc_y_scale($flag."_".$chrom,TRACK_HEIGHT);
    }
    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE)
				* $count) + TRACK_SPACE / 2;

    my $track_bottom = $track_baseline + 2*TRACK_HEIGHT;

    for my $point (@{$obj_top->{nonsliding}->{$chrom}}) {
	my ($pos,$val) = @{$point};

	my $left  = ($pos * $xscale) + TRACK_LEFT;
	my $top   = ($track_baseline + TRACK_HEIGHT) - ($val * $yscale_top);
	if ( $top < $track_baseline ) {
	    $top = $track_baseline;
	}
	if ($GD) {
	    $img->line($left,$top,$left,
		       $track_baseline + TRACK_HEIGHT,$fcolor);
	    } else {
		$img->line(x1=>$left, y1=>$top,
			   x2=>$left, y2=>$track_baseline + TRACK_HEIGHT,
			   id=>"$label-$pos-fwd",
			   stroke=>$fcolor,
			   fill=>$fcolor);
	    }
	# print STDERR join("\t",$left,$top),"\n";
    }

    for my $point (@{$obj_bot->{nonsliding}->{$chrom}}) {
	my ($pos,$val) = @{$point};

	my $left  = ($pos * $xscale) + TRACK_LEFT;
	my $top   = ($track_baseline + TRACK_HEIGHT) + ($val * $yscale_bot);
	if( $top > $track_bottom ) {
	    $top = $track_bottom;
	}
	if ($GD) {
	    $img->line($left,$top,$left,
		       $track_baseline + TRACK_HEIGHT,$rcolor);
	    } else {
		$img->line(x1=>$left, y1=>$top,
			   x2=>$left, y2=>$track_baseline + TRACK_HEIGHT,
			   id=>"$label-$pos-rev",
			   stroke=>$rcolor,
			   fill=>$rcolor);
	    }
	# print STDERR join("\t",$left,$top),"\n";
    }
    draw_yticks($track_baseline,
		$yscale_top,
		$yscale_bot,
		'5',$obj_top->{ylabel},
		0, # sig figures
		1, # two panel
		);

    draw_bounding($track_baseline,$obj_top,1);
    $CONFIG->{count} += 2;
    return;
}

sub plot_gmap {
  my ($gmap) = @_;
  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count};

  my $color = $settings->{black};
  my $black = $settings->{black};
  my $fontheight = $settings->{fontheight};
  $count--;

  my $min = $gmap->{$chrom}->{min};
  my $max = $gmap->{$chrom}->{max};
  my $range = $max + abs ($min);

  my $yscale = TRACK_HEIGHT / $range;
  my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;

  my @genes = @{$gmap->{$chrom}->{positions}};
  my @sorted = sort { $a->[1] <=> $b->[1] } @genes;
  my $features_seen;
  for my $gene (@sorted) {
    $features_seen++;
    my $gmap = $gene->[0] - $min;  # Set the mimumum gmap value to 0
    my $pmap = $gene->[1];
    my $left  = ($pmap * $xscale) + TRACK_LEFT;
    my $top   = $track_baseline + (TRACK_HEIGHT - ($gmap * $yscale));
    if ($GD) {
      $img->setPixel($left,$top,$color);
    } else {
      $img->circle(cx=>$left,cy=>$top,r=>'0.15',
		   id=>"gmap-$features_seen-$gmap-$pmap",
		   # id=>"$label-$bin",
		   stroke=>$color,
		   fill=>$color);
    }
    # print STDERR join("\t",$left,$top),"\n";
    # print STDERR $min,"\t",$max,"\n";
    $features_seen++;
  }
  # draw yticks
  my $width  = $CHROMS{$chrom}[1];
  my $tick_right  = TRACK_LEFT + ($xscale * $width);
  my $tick_left  = $tick_right - 4;

  if ($GD) {
    $img->line($tick_left,$track_baseline,$tick_right,$track_baseline,$settings->{black});
    $img->string(gdSmallFont,$tick_right+4,$track_baseline,,# - ($fontheight/2),
		 '25',$settings->{black});
    # The labels still might be shifted up inappropriately
    $img->line($tick_left,
	       $track_baseline + (TRACK_HEIGHT/2),$tick_right,
	       $track_baseline + (TRACK_HEIGHT/2),
	       $settings->{black});
    $img->string(gdSmallFont,$tick_right+4,$track_baseline + (TRACK_HEIGHT/2),# - ($fontheight/2),
		 '0',$settings->{black});
    $img->line($tick_left,
	       $track_baseline + (TRACK_HEIGHT),$tick_right,
	       $track_baseline + (TRACK_HEIGHT),
	       $settings->{black});
    $img->string(gdSmallFont,$tick_right+4,$track_baseline + (TRACK_HEIGHT),# - ($fontheight/2),
		 '-25',$settings->{black});
    # Draw the unit label on the right side of the chart
    $img->string(gdMediumBoldFont,$tick_right+14,$track_baseline + (TRACK_HEIGHT/2)
		 - ($fontheight/2),'cM',$settings->{black});
  } else {
    # 25 CM
    $img->line(x1=>$tick_left,y1=>$track_baseline,
	       x2=>$tick_right,y2=>$track_baseline,
	       id=>"gmap 25 CM ytick",
	       stroke => $settings->{black},
	       fill   => $settings->{black});
    $img->text(
	       style=> {
			'font' => 'Arial',
			'font-size' => LABEL_SIZE,
		       },
	       id=>"gmap 25 CM",
	       x=>$tick_right+4,
	       y=>$track_baseline)->cdata('25');
    # 0 CM
    $img->line(x1=>$tick_left,y1=>$track_baseline + (TRACK_HEIGHT/2),
	       x2=>$tick_right,y2=>$track_baseline + (TRACK_HEIGHT/2),
	       id=>"gmap 0 CM ytick",
	       stroke => $settings->{black},
	       fill   => $settings->{black});

    $img->text(
	       style=> {
			'font' => 'Arial',
			'font-size' => LABEL_SIZE,
		       },
	       id=>"gmap 0 CM",
	       x=>$tick_right+4,
	       y=>$track_baseline + (TRACK_HEIGHT/2)+2)->cdata('0');

    # -25 CM
    $img->line(x1=>$tick_left,y1=>$track_baseline + (TRACK_HEIGHT),
	       x2=>$tick_right,y2=>$track_baseline + (TRACK_HEIGHT),
	       id=>"gmap -25 CM ytick",
	       stroke => $settings->{black},
	       fill   => $settings->{black});

    $img->text(
	       style=> {
			'font' => 'Arial',
			'font-size' => LABEL_SIZE,
		       },
	       id=>"gmap -25 CM",
	       x=>$tick_right + 4,
	       y=>$track_baseline + (TRACK_HEIGHT)+2)->cdata('-25');

    # Draw the unit label on the right side of the chart
    $img->text(
	       style=> {
			'font' => 'Arial',
			'font-size' => LABEL_SIZE,
		       },
	       id=>"CM label",
	       x=>$tick_right + 14,
	       y=>$track_baseline + (TRACK_HEIGHT/2)+2)->cdata('cM');
  }
}

sub plot_posneg_numeric { 
    my ($label,$obj) = @_;
    my $PLOT_BY = 'average';
    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};

    my $color = $settings->{red};
    my $black = $settings->{red};

    # Calculate my own yscale
    my ($range);
    if ( $label eq 'tajimaD' ) {
	$range = 3;
    } else {
	$range = 1;
    }

    # The baseline for these plots will be 0, centered in the middle of the plot
    # with equidistant spacing on both sides...
    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;
    my $yscale = TRACK_HEIGHT / $range;

    # This is the old approach - plotting averages in bins across the chromosome
    my $already_plotted;
    if ( ! $already_plotted) {
	# This plots a seperate point for every ortholog pair
	for my $point (@{$obj->{nonsliding}->{$chrom}}) {
	    my ($pos,$val) = @{$point};

	    my $left  = ($pos * $xscale) + TRACK_LEFT;
	    my $top   = ($track_baseline + TRACK_HEIGHT) - ($val * $yscale);
	    if ( $top < $track_baseline ) {
		$top = $track_baseline;
	    }
	    if ($GD) {
		$img->line($left,$top,$left,$track_baseline + TRACK_HEIGHT,$color);
	    } else {
		$img->line(x1=>$left, y1=>$top,
			   x2=>$left, y2=>$track_baseline + TRACK_HEIGHT,
			   id=>"$label-$point",
			   stroke=>$color,
			   fill=>$color);
	    }
	    # print STDERR join("\t",$left,$top),"\n";
	}
    }
    # Draw a grey line at 1.0 all the way across the graph
    # for Tajima D

    my $width  = $CHROMS{$chrom}[1];
    my $right  = TRACK_LEFT + ($xscale * $width) + 2;
    my $location = ($track_baseline + TRACK_HEIGHT);
    if ($GD) {
	$img->line(TRACK_LEFT,$location,$right,$location,$settings->{grey});
    } else {
	$img->line(x1=>TRACK_LEFT,
		   y1=>$location,
		   x2=>$right,
		   y2=>$location,
		   id=>"$label-significance divider",
		   stroke => $settings->{grey},
		   fill   => $settings->{grey});
    }

    draw_yticks($track_baseline,$yscale,$yscale,'5',$obj->{ylabel}, 
		1, # sig figs are 1
		1, # two-paned pos/neg plot
		);
    draw_bounding($track_baseline,$obj,1); # the third parameter indicates 
                                           # this is a pos and neg pane
    $CONFIG->{count}++;
    return;
}

sub plot_numeric { 
    my ($label,$obj,$range,$color,$back) = @_;
    my $PLOT_BY = 'average';
    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};
    
    $color = $color ? $settings->{$color} : $settings->{red};
    $back = $back ? $settings->{$back} : $settings->{black};
    
    # Calculate my own yscale
    unless( defined $range ) {
	if ($label eq 'pi') {
	    $range = 0.06;
	} elsif ($label eq 'segsites') {
	    $range = 100;
	} else {
	    $range = 1;
	}
    }

    # The baseline for these plots will be 0, centered in the middle of the plot
    # with equidistant spacing on both sides...
    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;
    my $yscale = TRACK_HEIGHT / $range;
    
    # This is the old approach - plotting averages in bins across the chromosome
    my $already_plotted;
    if ( ! $already_plotted) {
	# This plots a seperate point for every ortholog pair
	for my $point (@{$obj->{nonsliding}->{$chrom}}) {
	    my ($pos,$val) = @{$point};

	    my $left  = ($pos * $xscale) + TRACK_LEFT;
	    my $top   = ($track_baseline + TRACK_HEIGHT) - ($val * $yscale);
	    if ( $top < $track_baseline ) {
		$top = $track_baseline;
	    }
	    if ($GD) {
		$img->line($left,$top,$left,$track_baseline + TRACK_HEIGHT,$color);
	    } else {
		$img->line(x1=>$left, y1=>$top,
			   x2=>$left, y2=>$track_baseline + TRACK_HEIGHT,
			   id=>"$label-$point",
			   stroke=>$color,
			   fill=>$color);
	    }
	    # print STDERR join("\t",$left,$top),"\n";
	}
    }

    draw_yticks($track_baseline,$yscale,$yscale,'5',
		$obj->{ylabel}, 
		$label =~ /^pi\./ ? 2 : 1, # sig figs calculation
		0, # single pane
		);
    draw_bounding($track_baseline,$obj);
    $CONFIG->{count}++;
    return;
}

sub plot_blocks {
    my ($obj,$icolor,$ricolor) = @_;
    my $fcolor = $settings->{$icolor};
    my $rcolor = $settings->{$ricolor};
    my $black = $settings->{black};

    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};

    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;
    if( ! exists $obj->{$chrom} ) {
	warn( " no chrom $chrom available\n");
    }
    my @blocks = @{$obj->{$chrom} || []};

    # assume that the strand is the majority of blocks 

    my %majstrand;
    { 
	my %c;  
	for my $block (@blocks) {	  
	    my ($start,$end,$target,$tlength,$tstrand) = @$block;
	    $c{$target}->{$tstrand} += $tlength;
	}
	while( my ($target,$pieces) = each %c ) {
	    ($majstrand{$target}) = sort { $pieces->{$b} <=> $pieces->{$a} } keys %{$pieces};
	}
    }

    my $total;
    for my $block (@blocks) {
	$total++;
	my ($start,$stop,$target,$tlength,$tstrand) = @$block;
	($start,$stop) = ($stop,$start) if ($start > $stop);
	my $left   = ($start * $xscale) + TRACK_LEFT;
	my $right  = ($stop * $xscale) + TRACK_LEFT;
	my $top    = $track_baseline;
	my $bottom = $track_baseline + TRACK_HEIGHT;
	if ($GD) {
	    $img->filledRectangle($left,$top+1,$right-1,$bottom-1,
				  $tstrand eq $majstrand{$target} ? 
				  $fcolor : $rcolor);
	} else {
	    $img->rectangle(x=>$left,y=>$top,
			    width  =>$right-1-$left,
			    height =>$bottom-$top,
			    id     =>$obj->{dumped_file}."-$total-" . $start .'-' . $stop,
			    stroke => ($tstrand eq $majstrand{$target} ? 
				       $fcolor : $rcolor),
			    fill   => ($tstrand eq $majstrand{$target} ? 
				       $fcolor : $rcolor));
	}	 
	# print STDERR join("\t",$left-1,$top,$right+1,$bottom),"\n";
    }

    draw_bounding($track_baseline,$obj);
    $CONFIG->{count}++;
    return;
}

sub plot_fish_blocks {
    my ($obj1,$obj2,$icolor1,$icolor2) = @_;
    my $color1 = $settings->{$icolor1};
    my $color2 = $settings->{$icolor2};
    my $black = $settings->{black};

    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};

    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;
    if( ! exists $obj1->{$chrom} ) {
	warn( " no chrom $chrom available\n");
    }
    my @blocks = @{$obj1->{$chrom} || []};

    # assume that the strand is the majority of blocks 

    my %majstrand;
    { 
	my %c;  
	for my $block (@blocks) {	  
	    my ($start,$end,$target,$tlength,$tstrand) = @$block;
	    $c{$target}->{$tstrand} += $tlength;
	}
	while( my ($target,$pieces) = each %c ) {
	    ($majstrand{$target}) = sort { $pieces->{$b} <=> $pieces->{$a} } keys %{$pieces};
	}
    }

    my $total;
    for my $block (@blocks) {
	$total++;
	my ($start,$stop,$target,$tlength,$tstrand) = @$block;
	($start,$stop) = ($stop,$start) if ($start > $stop);
	my $left   = ($start * $xscale) + TRACK_LEFT;
	my $right  = ($stop * $xscale) + TRACK_LEFT;
	my $top    = $track_baseline;
	my $bottom = $track_baseline + TRACK_HEIGHT;
	if ($GD) {
	    $img->filledRectangle($left,$top+1,$right-1,$bottom-1,$color1);
	} else {
	    $img->rectangle(x=>$left,y=>$top,
			    width  =>$right-1-$left,
			    height =>$bottom-$top,
			    id     =>$obj1->{dumped_file}."-$total-" . $start .'-' . $stop,
			    stroke => $color1,
			    fill   => $color1);
	}	 
	# print STDERR join("\t",$left-1,$top,$right+1,$bottom),"\n";
    }
    @blocks = @{$obj2->{$chrom} || []};

    for my $block (@blocks) {
	$total++;
	my ($start,$stop,$target,$tlength,$tstrand) = @$block;
	($start,$stop) = ($stop,$start) if ($start > $stop);
	my $left   = ($start * $xscale) + TRACK_LEFT;
	my $right  = ($stop * $xscale) + TRACK_LEFT;
	my $top    = $track_baseline;
	my $bottom = $track_baseline + TRACK_HEIGHT;
	if ($GD) {
	    $img->filledRectangle($left,$top+1,$right-1,$bottom-1,$color2);
	} else {
	    $img->rectangle(x=>$left,y=>$top,
			    width  =>$right-1-$left,
			    height =>$bottom-$top,
			    id     =>$obj2->{dumped_file}."-$total-" . $start .'-' . $stop,
			    stroke => $color2,
			    fill   => $color2);
	}	 
	# print STDERR join("\t",$left-1,$top,$right+1,$bottom),"\n";
    }

    draw_bounding($track_baseline,$obj2);
    $CONFIG->{count}++;
    return;
}

sub plot_centromere_telomere {
    my ($obj,$icencolor,$itelcolor) = @_;
    my $cencolor = $settings->{$icencolor};
    my $telcolor = $settings->{$itelcolor};

    my $black = $settings->{black};

    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};

    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;
    if( ! exists $obj->{$chrom} ) {
	warn( " no chrom $chrom available\n");
    }
    my @blocks = @{$obj->{$chrom} || []};

    my $total;
    for my $block (@blocks) {
	$total++;
	my ($start,$stop,$type) = @$block;
	my $color_stroke = $type =~ /^CEN/i ? $cencolor : $telcolor;
	($start,$stop) = ($stop,$start) if ($start > $stop);
	my $left   = ($start * $xscale) + TRACK_LEFT;
	my $right  = ($stop * $xscale) + TRACK_LEFT;
	my $top    = $track_baseline;
	my $bottom = $track_baseline + TRACK_HEIGHT / 2;
	if ($GD) {
	    $img->filledRectangle($left,$top+1,$right-1,$bottom-1,
				  $color_stroke);
	} else {
	    $img->rectangle(x=>$left,y=>$top,
			    width  => $right-1-$left,
			    height => $bottom-$top,
			    id     => $type."-$total-" . $start .'-' . $stop,
			    stroke => $color_stroke,
			    fill   => $color_stroke);	    
	}	 
	# print STDERR join("\t",$left-1,$top,$right+1,$bottom),"\n";
    }
    draw_bounding($track_baseline,$obj,undef,1);
    $CONFIG->{count} += 0.75;
    return;
}

sub plot_recombination_rates {
    my ($obj,$hicolor,$cicolor,$neutralcolor,$markers) = @_;
    
    my %colormap = ( 'HOT' => $settings->{$hicolor},
		     'COLD' => $settings->{$cicolor},
		     'NEUTRAL' => $settings->{$neutralcolor},
		     'UNSCORED' => $settings->{'white'});
    my $black = $settings->{black};

    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};

    my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;
    if( ! exists $obj->{$chrom} ) {
	warn( " no chrom $chrom available\n");
    }
    my @blocks = @{$obj->{$chrom} || []};
    my @ssrs;
    if( defined $markers ) {
	@ssrs = @{$markers->{$chrom} || []};
    }
    # assume that the strand is the majority of blocks 

    my $total;
    for my $block (@blocks) {
	$total++;
	my ($start,$stop,$target,$tlength,$hot_cold) = @$block;
	# hot_cold will be -1,0,1
	($start,$stop) = ($stop,$start) if ($start > $stop);
	my $left   = ($start * $xscale) + TRACK_LEFT;
	my $right  = ($stop * $xscale) + TRACK_LEFT;
	my $top    = $track_baseline;
	my $bottom = $track_baseline + TRACK_HEIGHT;
	if ($GD) {
	    $img->filledRectangle($left,$top+1,$right-1,$bottom-1,
				  $colormap{$hot_cold});
	} else {
	    $img->rectangle(x=>$left,y=>$top,
			    width  => $right-1-$left,
			    height => $bottom-$top,
			    id     => $target."-$total-" . $start .'-' . $stop,
			    stroke => $colormap{$hot_cold},
			    fill   => $colormap{$hot_cold});
	}	 
	# print STDERR join("\t",$left-1,$top,$right+1,$bottom),"\n";
    }
    $total = 0;
    for my $marker ( @ssrs ) {
	$total++;
	my ($start,$stop) = @$marker;
	($start,$stop) = ($stop,$start) if ($start > $stop);
	my $left   = ($start * $xscale) + TRACK_LEFT;
#	my $right  = ($stop * $xscale) + TRACK_LEFT;
	my $right  = $left + 1;
#	warn("left is $left, right is $right\n for $start..$stop\n");
	my $top    = $track_baseline;
	my $bottom = $track_baseline + TRACK_HEIGHT;
	if ($GD) {
	    $img->filledRectangle($left,$top+1,$right-1,$bottom-1,
				  $black);
	} else {
	    $img->rectangle(x=>$left,y=>$top,
			    width  =>$right - 0.90 - $left,
			    height =>$bottom-$top,
			    id     =>"SSR-$total-" . $start .'-' . $stop,
			    stroke => $black,
			    fill   => $black);
	}	 
    }
    draw_bounding($track_baseline,$obj);
    $CONFIG->{count}++;
    return;
}


# Draw a bounding box for this track.
sub draw_bounding {
  my ($top,$obj,$two_pane,$halfpane) = @_;
  my $track_label = $obj->{label};
  my $units       = $obj->{units};
  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count}; 

  my $width  = $CHROMS{$chrom}[1];
  my $left   = TRACK_LEFT;
  my $right  = TRACK_LEFT + ($xscale * $width);
  
  my $bottom = $top + TRACK_HEIGHT;
  if( $halfpane ) {
      $bottom -= TRACK_HEIGHT / 2;
  } elsif( $two_pane ) {
      $bottom += TRACK_HEIGHT;
  }

  if ($GD) {
    $img->rectangle($left,$top,$right,$bottom,$settings->{black});
  } else {
    $img->rectangle(x     =>$left,
		    y     =>$top,
		    width =>$right-$left,
		    height=>$bottom-$top,
		    id  => $track_label,
		    stroke=>$settings->{black},
		    'fill-opacity' => 0);
  }

  # Draw xscale tickmarks every million basepairs
  my $tick_top     = $top + TRACK_HEIGHT - 4;
  my $tick_bottom  = $top + TRACK_HEIGHT;

  if( $two_pane || $halfpane ) { 
      $tick_bottom  = $bottom;
      $tick_top     = $bottom - 4;
  } 
  my $chrom_length = $CHROMS{$chrom}[1];
  my $i;
  for (my $i=0;$i<=$chrom_length;$i+= TICK_SCALE ){
    my $left = TRACK_LEFT + ($i * $xscale);
    my $label = $i / TICK_SCALE;

    if ($GD) {
      $img->line($left,$tick_top,$left,$tick_bottom,$settings->{black});
      $img->string(gdSmallFont,
		   $left-length($label),
		   $tick_bottom + 2,$label,$settings->{black});
    } else {
      $img->line(x1=>$left,y1=>$tick_top,
		 x2=>$left,y2=>$tick_bottom,
		 id=>"$track_label-$label MBp tick",
		 stroke => $settings->{black},
		 fill   => $settings->{black});

      $img->text(
		 style=> {
			  'font' => 'Arial',
			  'font-size' => LABEL_SIZE,
			  'color'=> $settings->{black},
			 },
		 id=>"$track_label-$label MBp",
		 x=>$left-3,
		 y=>$tick_bottom + 18)->cdata($label);
    }
  }
  # Generate the track label
  # On the right side of the box...
  #  my $label_left = 
  #    $left
  #      + ($xscale * $chrom_length) - (length($track_label) * $fontwidth) - 10;
  # or on the left side of the box...
  # my $label_left = TRACK_LEFT + 4;
  # my $track_top = $top;

  # ...or just above the box
  my $label_left = TRACK_LEFT + 4;
  my $label_top = $top - (gdMediumBoldFont->height) - 4;

  if ($GD) {
    $img->string(gdMediumBoldFont,
		 $label_left,
		 $label_top + 3,$track_label,$settings->{black});
  } else {
    $img->text(
	       style=> {
			'font' => 'Arial',
			'font-size' => 22,
		       },
	       id=>"$track_label-main label",
	       x=>$label_left,
	       y=>$top - 8)->cdata($track_label);
  }

  # ...or on the left side of the box
  # Center everything based on the longest label
  #my @labels = sort { length $b <=> length $a } values %ylabels;
  #my $longest = (length $labels[0]) * gdMediumBoldFont->width;

  #my $label_left = 4 + (($longest - (length $track_label) * gdMediumBoldFont->width)/2);
  #my $label_top = $bottom - (TRACK_HEIGHT/2) - (gdMediumBoldFont->height);

  #$img->string(gdMediumBoldFont,
  # 	       $label_left,
  #       $label_top + 3,$track_label,$settings->{black});
}



# This is a factory for new objects...
sub establish_image {
  if ($GD) {
    my $img = GD::Image->new(WIDTH,
			     TOTAL_TRACKS * (TRACK_HEIGHT + TRACK_SPACE));
    return $img;
  } else {
    my $img = SVG->new('width' => WIDTH,
		       'height'=> TOTAL_TRACKS * (TRACK_HEIGHT + TRACK_SPACE) + 25);
    return $img;
  }
}

sub establish_settings {
  my $img = shift;
  my $settings = {};
  if ($GD) {
      # Establish some colors and dimensional values.
      $settings = { white  => $img->colorAllocate(255,255,255),
		    black  => $img->colorAllocate(0,0,0),
		    purple => $img->colorAllocate(75,000,130),
		    grey   => $img->colorAllocate(230,230,230),
		    pink   => $img->colorAllocate(204,000,204),
		    sea    => $img->colorAllocate(000,102,102),
		    orange => $img->colorAllocate(255,153,000),
		    darkorange => $img->colorAllocate(255,102,000),
		    darkblue => $img->colorAllocate(6,38,111),
		    medblue  => $img->colorAllocate(42,68,128),
		    darkgreen => $img->colorAllocate(0,51,0),
		    lightgreen => $img->colorAllocate(0,153,0),
		    darkred => $img->colorAllocate(153,0,0),
		    lightred => $img->colorAllocate(204,51,51),
		    seablue => $img->colorAllocate(0,102,204),
		    limegreen => $img->colorAllocate(51,255,51),
		    darkgrey => $img->colorAllocate(170,170,170),
		    # *** Primary Color:

		    # var. 1 = #133AAC = rgb(19,58,172)
		    # var. 2 = #2B4181 = rgb(43,65,129)
		    # var. 3 = #062170 = rgb(6,33,112)
		    # var. 4 = #476BD6 = rgb(71,107,214)
		    # var. 5 = #6D87D6 = rgb(109,135,214)
		    varP1 => $img->colorAllocate(19,58,172),
		    varP2 => $img->colorAllocate(43,65,129),
		    varP3 => $img->colorAllocate(6,33,112),
		    varP4 => $img->colorAllocate(71,107,214),
		    varP5 => $img->colorAllocate(109,135,214),
		    # *** Secondary Color A:

		    # var. 1 = #FFE100 = rgb(255,225,0)
		    # var. 2 = #BFAE30 = rgb(191,174,48)
		    # var. 3 = #A69200 = rgb(166,146,0)
		    # var. 4 = #FFE840 = rgb(255,232,64)
		    # var. 5 = #FFEE73 = rgb(255,238,115)

		    varA1 => $img->colorAllocate(255,255,0),
		    varA2 => $img->colorAllocate(191,174,48),
		    varA3 => $img->colorAllocate(166,146,0),
		    varA4 => $img->colorAllocate(255,232,64),
		    varA5 => $img->colorAllocate(255,238,115),
		    # *** Secondary Color B:

		    # var. 1 = #FF6A00 = rgb(255,106,0)
		    # var. 2 = #BF6B30 = rgb(191,107,48)
		    # var. 3 = #A64500 = rgb(166,69,0)
		    # var. 4 = #FF8F40 = rgb(255,143,64)
		    # var. 5 = #FFAD73 = rgb(255,173,115)		    

		    varB1 => $img->colorAllocate(255,106,0),
		    varB2 => $img->colorAllocate(191,107,48),
		    varB3 => $img->colorAllocate(166,69,0),
		    varB4 => $img->colorAllocate(255,143,64),
		    varB5 => $img->colorAllocate(255,173,115),

		    # Good colors
		    red    => $img->colorAllocate(255,0,0),
		    blue   => $img->colorAllocate(000,000,255),
		    green  => $img->colorAllocate(000,255,000),
		    aqua   => $img->colorAllocate(000,255,204),
		    rosybrown => $img->colorAllocate(188,143,143),
		    indianred=> $img->colorAllocate(205,92,92),
		    brown  => $img->colorAllocate(139,69,19),
		    olive  => $img->colorAllocate(107,142,35),
		    steelblue => $img->colorAllocate(70,130,180),
		    goldenrod => $img->colorAllocate(238,221,130),
		    slateblue => $img->colorAllocate(106,90,205),
		    forestgreen=> $img->colorAllocate(34,139,34),

		    fontwidth  => gdMediumBoldFont->width,
		    fontheight => gdMediumBoldFont->height,
		};
  } else {
    $settings = { white  => 'rgb(255,255,255)',
		  black  => 'rgb(0,0,0)',
		  purple => 'rgb(200,100,200)',
		  grey   => 'rgb(230,230,230)',
		  darkgrey => 'rgb(170,170,170)',
		  pink   => 'rgb(204,000,204)',
		  sea    => 'rgb(000,102,102)',
		  orange => 'rgb(255,153,000)',
		  darkorange => 'rgb(255,102,000)',
		  limegreen => 'rgb(51,255,51)',
		  lightgreen=> 'rgb(0,153,0)',
		  lightred => 'rgb(204,51,51)',
		    
		    # *** Primary Color:

		    # var. 1 = #133AAC = rgb(19,58,172)
		    # var. 2 = #2B4181 = rgb(43,65,129)
		    # var. 3 = #062170 = rgb(6,33,112)
		    # var. 4 = #476BD6 = rgb(71,107,214)
		    # var. 5 = #6D87D6 = rgb(109,135,214)
		    varP1 => 'rgb(19,58,172)',
		    varP2 => 'rgb(43,65,129)',
		    varP3 => 'rgb(6,33,112)',
		    varP4 => 'rgb(71,107,214)',
		    varP5 => 'rgb(109,135,214)',
		    # *** Secondary Color A:

		    # var. 1 = #FFE100 = rgb(255,225,0)
		    # var. 2 = #BFAE30 = rgb(191,174,48)
		    # var. 3 = #A69200 = rgb(166,146,0)
		    # var. 4 = #FFE840 = rgb(255,232,64)
		    # var. 5 = #FFEE73 = rgb(255,238,115)

		    varA1 => 'rgb(255,255,0)',
		    varA2 => 'rgb(191,174,48)',
		    varA3 => 'rgb(166,146,0)',
		    varA4 => 'rgb(255,232,64)',
		    varA5 => 'rgb(255,238,115)',
		    # *** Secondary Color B:

		    # var. 1 = #FF6A00 = rgb(255,106,0)
		    # var. 2 = #BF6B30 = rgb(191,107,48)
		    # var. 3 = #A64500 = rgb(166,69,0)
		    # var. 4 = #FF8F40 = rgb(255,143,64)
		  # var. 5 = #FFAD73 = rgb(255,173,115)		    
		  
		  varB1 => 'rgb(255,106,0)',
		  varB2 => 'rgb(191,107,48)',
		  varB3 => 'rgb(166,69,0)',
		  varB4 => 'rgb(255,143,64)',
		  varB5 => 'rgb(255,173,115)',
		  # Good colors
		  red    => 'rgb(255,0,0)',
		  blue   => 'rgb(000,000,255)',
		  green  => 'rgb(000,255,000)',
		  aqua   => 'rgb(000,255,204)',
		  brown  => 'rgb(139,69,19)',
		  olive  => 'rgb(107,142,35)',
		  rosybrown   => 'rgb(188,143,143)',
		  steelblue   => 'rgb(70,130,180)',
		  goldenrod   => 'rgb(238,221,130)',
		  slateblue   => 'rgb(106,90,205)',
		  forestgreen => 'rgb(34,139,34)',
		  darkblue => 'rgb(6,38,111)',
		  medblue  => 'rgb(42,68,128)',
		  darkgreen => 'rgb(0,51,0)',
		  darkred => 'rgb(153,0,0)',
		  seablue  => 'rgb(0,102,204)',

		};
  }
  return $settings;
}


package Windows;

use List::Util qw(sum max);
# constants for sliding windows
use constant WINDOW        => 50_000;
use constant STEP          => 1_000;

sub new {
  my ($self,$label) = @_;
  my $this = bless {},$self;
  $this->{dumped_file} = $label;
  $this->{label}  = $labels{$label};
  $this->{ylabel} = $ylabels{$label};
  $this->{units}  = $units{$label};
  $this->{fh} = undef;
  return $this;
}


sub read_dumped_file {
  my $self = shift;
  my $file = $self->{dumped_file};
  open(my $fh => $file) || die "$!\n";
  while (<$fh>) {
    chomp;
    # Commented lines are to be loaded directly into the object
    # I'm really just trying to recover the max_y values here...
    if (/^\#/) {
      # parse out the field and the value
      $_ =~ /^\#(.*)=(.*)/;
      my $field = $1;
      my $val   = $2;
      $self->{$field} = $val;
    }
    my ($chrom,$bin,$val) = split("\t");
    # Redundancy in the data structure since I don't know what was printed out
    $self->{$chrom}->{$bin}->{total} = $val;
    $self->{$chrom}->{$bin}->{average} = $val;
  }
  $self->{fh} = $fh;
  1;
}


sub parse {
    my ($self,$file,$group_by,$bin_by,$save_by, $filter) = @_;
    my $cols = fetch_columns($file);
    for my $col ( $group_by, $bin_by, $save_by) {
	if( defined $col && ! exists $cols->{$col}) {
	    die("cannot find column $col in $file\n");
	}
    }
    warn( "parsing: $file ...\n");
    my $positions = {};
    if ($USE_CACHED) {
	$self->read_dumped_file();
	return;
    } else {
	open($self->{fh} => $file) or die "$! $file\n";;
    }
    my $fh = $self->{fh};
    while (<$fh>) {
	chomp;
	# Skip comments
	next if (/^\#/);
	# Fetch the position of the bin_by and group_by
	# columns in the fields array
	my @fields    = split("\t",$_);
	my $bin_val   = $fields[$cols->{$bin_by}];
	my $group_val = $fields[$cols->{$group_by}];

	# Is the value out of range?
	if ($bin_val > $CHROMS{$group_val}[1]) {	
	    print STDERR "Positional value out of range for chromosome $group_val...$bin_val\t",
	    $CHROMS{$group_val}[1],"\n";
	    next;
	}
	# Save both the physical position to map and the value...
	# (these may or may not be the same thing!)
	# This is normally used for things like KaKs values, 
	# where I am plotting values (and not just sums of occurences) 
	# against the physical position
	my $value = 0;
	if( defined $save_by ) {
	    $value  = eval { $fields[$cols->{$save_by}] };
	    if( ! defined $filter || &$filter($value))  {
		push (@{$self->{nonsliding}->{$group_val}},[$bin_val,$value]);
	    } else {
		next;
	    }
	}  else {
	    $value = 1;
	    # Okay, I didn't find a save_by value. Just plotting by the bin_value
	}
	$value   ||= $bin_val;
	push (@{$positions->{$group_val}},[$bin_val,$value]);
    }
    $self->sliding_windows($positions);
    return;
}

sub parse_blocks {
    my ($self,$file,$chrom_col,$start_col,$end_col,$target_col) = @_;
    print STDERR "parsing: $file ...\n";
    my $cols = fetch_columns($file);
    for my $col ( $chrom_col,$start_col,$end_col,$target_col) { 
	if( defined $col && ! exists $cols->{$col}) {
	    die("cannot find column $col in $file\n");
	}
    }
    my $positions = {};
    open IN,$file or die "$! $file\n";
    while (<IN>) {
	chomp;
	# Skip comments
	next if (/^\#/);
	my @fields = split(/\t/,$_);
	my $chrom = $fields[$cols->{$chrom_col}];
	my $start = $fields[$cols->{$start_col}];
	my $end   = $fields[$cols->{$end_col}];
	my $target = $fields[$cols->{$target_col}];
	push (@{$self->{$chrom}},[$start,$end,$target,abs($end-$start),'+']);
    }
    return;
}

sub parse_mercator_blocks {
  my ($self,$file) = @_;
  print STDERR "parsing: $file ...\n";

  my $positions = {};
  open IN,$file or die "$! $file\n";
  while (<IN>) {
    chomp;
    # Skip comments
    next if (/^\#/);
    my($chrom,$start,$end,$strand,$target,$tstart,$tend,$tstrand) = split "\t";
    # Fetch the position of the bin_by and group_by
    # columns in the fields array
    push (@{$self->{$chrom}},[$start,$end,$target,abs($tend-$tstart),$tstrand]);
  }
  return;
}

sub parse_range {
  my ($self,$file,$group_by,$start,$end) = @_;
  print STDERR "parsing: $file ...\n";
  my $cols = fetch_columns($file);
  for my $col ( $group_by, $start, $end) {
      if( defined $col && ! exists $cols->{$col}) {
	  die("cannot find column $col in $file\n");
      }
  }
  my $positions = {};
  open IN,$file or die "$! $file\n";
  while (<IN>) {
      chomp;
      # Skip comments
      next if (/^\#/);
      # Fetch the position of the bin_by and group_by
      # columns in the fields array
      my @fields    = split("\t",$_);
      my $group_val = $fields[$cols->{$group_by}];
      my $start_val = $fields[$cols->{$start}];
      my $end_val   = $fields[$cols->{$end}];

      # Fetch the position of the bin_by and group_by
      # columns in the fields array
      push @{$self->{$group_val}},[$start_val,$end_val];
  }
  return;
}

sub parse_recombination_rates {
  my ($self,$file) = @_;
  print STDERR "parsing: $file ...\n";
  my $cols = fetch_columns($file);
  for my $col ( qw(start stop chrom type designation) ) {
      if( defined $col && ! exists $cols->{$col}) {
	  die("cannot find column $col in $file\n");
      }
  }  
  my $positions = {};
  open IN,$file or die "$! $file\n";
  while (<IN>) {
    chomp;
    # Skip comments
    next if (/^\#/);
    my @fields    = split("\t",$_);
    my $start   = $fields[$cols->{'start'}];
    my $stop    = $fields[$cols->{'stop'}];
    my $chrom   = $fields[$cols->{'chrom'}];
    my $target  = $fields[$cols->{'type'}];
    my $desig   = $fields[$cols->{'designation'}];
    push (@{$self->{$chrom}},[$start,$stop,$target,abs($stop-$start),
			      $desig]);
  }
  return;
}

sub parse_gmap {
  my ($self,$file) = @_;
  open IN,$file or die "$! $file\n";
  my ($min,$max);
  while (<IN>) {
    chomp;
    # Skip comments
    next if (/^\#/);
    # Fetch the position of the bin_by and group_by
    # columns in the fields array
    # Strucutre is gmap position,pmap position for each gene
    my @fields    = split("\t");
    next if ($fields[3] eq '');

    push (@{$self->{$fields[0]}->{positions}},[$fields[3],$fields[4]]);
    $min = ($fields[3] < $min) ? $fields[3] : $min;
    $max = ($fields[3] > $max) ? $fields[3] : $max;
    $self->{$fields[0]}->{min} = $min;
    $self->{$fields[0]}->{max} = $max;
  }
  return;
}

# Which windows does this fall into?
sub sliding_windows {
    my ($self,$positions) = @_;

    print STDERR "   binning...\n";
    for my $group_value (keys %$positions) {
	# Get the maximal limit (ie the highest scoring feature)
	my @positions = sort { $a->[0] <=> $b->[0] } @{$positions->{$group_value}};
	for my $temp (@positions) {
	    my ($fstart,$value) = @$temp;
	    $self->stuff($fstart,$value,$positions[-1]->[0],$group_value);
	}
    }
    # calculate the average and total values in each bin
    $self->calc_averages();
}




# This is a new attempt at a vastly more efficient 
# calculated sliding window algorithm
sub stuff {
    my ($self,$fstart,$value,$max,$key) = @_;

    # Calculate the center bin position based on the STEP size,
    # not the WINDOW size...
    my $center_bin_index = int ($fstart / STEP);

    # How many steps are there per window?
    my $steps = WINDOW / STEP;

    # One approach: the feature should fall into equal bins on both sides
    # This centers each window on the bin point - maybe not quite accurate...
    # my $start = $center_bin_index - ($steps/2);
    # my $stop  = $center_bin_index + ($steps/2);

    # Second approach: the center_bin_index is really the LAST bin
    # that should contain the feature (that is, it's the upper limit
    # for bins containing the feature).
    my $start = $center_bin_index - $steps; # 0-based indexing
    my $stop  = $center_bin_index + 1;

    # Is $start less than 0?
    $start = ($start < 0) ? 0 : $start;

    # Non-overlapping bins...
    if ($steps == 1) {
	$start = $center_bin_index;
	$stop  = $center_bin_index;
    }
    for (my $i=$start;$i<=$stop;$i+=1) {
	my $bin = $i * STEP;
	push (@{$self->{groups}->{$key}->{$bin}->{values}},$value);
#	print STDERR "\t",$i,"\t",$bin,"\t",$value,,"\t",$fstart,"\n";
    }
    return;
}

# Calculate an average for each bin or a total value.
sub calc_averages {
    my $self = shift;
    warn("   calculating average...\n");
    my ($max_average,$max_total,$max_value_all) = (0,0,0);
    
    # maybe should also calculate 1SD to help in setting an upper limit cutoff?
    for my $group_by ($self->fetch_groups()) {	
	warn("group is empty") unless defined $group_by;
	my $bins = $self->fetch_bins($group_by);
	my $max_value = 0;
	for my $bin (@$bins) {
	    my $all = $self->fetch_values($group_by,$bin);

	    my $sum = &sum (@$all);
	    my $avg = $sum / scalar @$all;
	    $self->{groups}->{$group_by}->{$bin}->{average} = $avg;
	    $self->{groups}->{$group_by}->{$bin}->{total} = scalar @$all;
	    
	    $max_value = max(@$all, $max_value);
	    $max_total   = max( scalar @$all, $max_total);
	    $max_average = max($avg,$max_average);
	}
	$self->{'max_y_value_'.$group_by}   = $max_value;
	$max_value_all = max($max_value, $max_value_all);
    }
    $self->{'max_y_total'} = $max_total || 1;
    $self->{'max_y_value'} = $max_value_all || 1;
    $self->{'max_y_average'} = $max_average || 1;
}

sub fetch_columns {
  my $file = shift;
  open(IN,$file)|| die "$! $file\n";
  my %cols;
  while (<IN>) {
    chomp;
    my $line = $_;
    $line =~ s/\#//g;
    my @cols = split("\t",$line);
    my $pos = 0;
    for (@cols) {
      $cols{$_} = $pos;
      $pos++;
    }
    last;
  }
  return \%cols;
}


# Normalization
sub normalize {
    my ($numerator,$denominator) = @_;

    # Iterate through all the bins in the numerator
    # looking up the corresponding values in the denominator.

    my $max;
    print STDERR "   normalizing...\n";
    for my $group ($numerator->fetch_groups()) {
	my $bins = $numerator->fetch_bins($group);
	for my $bin (@$bins) {
	    my $total = $numerator->total($group,$bin);
	    my $denom_total = $denominator->total($group,$bin);
	    my $normalized = eval { $total / $denom_total };
	    $normalized ||= '0';
	    $max ||= $normalized;
	    $max = max($normalized,$max);
	    $numerator->{groups}->{$group}->{$bin}->{normalized} = $normalized;
	}
    }
    $numerator->{max_y_normalized} = $max;
}

# Data access methods
sub fetch_groups {  return keys %{shift->{groups}}; }

sub fetch_bins {
  my ($self,$group) = @_;
  my @bins = keys %{$self->{groups}->{$group}};
  my @sorted = sort { $a <=> $b } @bins;
  return \@sorted;
}
sub fetch_values {
  my ($self,$group,$bin) = @_;
  my @vals = @{$self->{groups}->{$group}->{$bin}->{values}};
  return \@vals;
}

sub total   {  
  my ($self,$group,$bin) = @_;
  return $self->{groups}->{$group}->{$bin}->{total};
}

sub average {  
  my ($self,$group,$bin) = @_;
  return $self->{groups}->{$group}->{$bin}->{average};
}

sub normalized {  
  my ($self,$group,$bin) = @_;
  return $self->{groups}->{$group}->{$bin}->{normalized};
}



sub print_info {
  my ($self,$tag,$field) = @_;
  return unless ($USE_CACHED);
  $field ||= 'total';
  my $file = $self->{dumped_file};
  open OUT,">dumped_values/$file.out";
  for my $group ($self->fetch_groups()) {

    # Print out some header information
    print OUT "#chromosome=$group\n";
    print OUT "#field_content=$field\n";
    print OUT "#max_y_average=" . $self->{max_y_average} . "\n";
    print OUT "#max_y_total=" . $self->{max_y_total} . "\n";
    print OUT "#max_y_normalized=" . $self->{max_y_normalized} . "\n";
    my $bins = $self->fetch_bins($group);
    for my $bin (@$bins) {
      my $normalized = $self->normalized($group,$bin);
      my $total = $self->$field($group,$bin);
      print OUT $group,"\t",$bin,"\t",$total,"\n";
    }
  }
  close OUT;
}


sub calc_y_scale {
  my ($self,$tag,$height) = @_;
  my $max = $self->{'max_y_' . $tag};
  unless( defined $max ) { warn("no max for ".'max_y_' . $tag."\n") }
  return 1 unless $max;
  return $height / $max;
}



1;
