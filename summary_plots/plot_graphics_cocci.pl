#!/usr/bin/perl -w
# $Id: plot_graphics.pl 250 2008-12-12 22:20:20Z stajich $

use strict;
use File::Spec;
use Bio::SeqIO;
use List::Util qw(sum max);
use GD qw(gdGiantFont gdLargeFont gdMediumBoldFont gdSmallFont gdTinyFont);
use SVG;

# SHOULD WE GENERATE THE PLOT WITH GD OR SVG?
# This should just be a command line parameter
my $GD  = 1;
my $SVG = 0;

use constant TICK_SCALE => 1_000_000; # draw tickmark every 1M bp

# recalculate values from all files...
my $USE_CACHED = 0;

my $DIR = '/projects/coccidiodies/variation/plot';

#celegans data
# DATA SOURCES
use constant GENES   => 'ci_gene_summary.tab';
#use constant GMAP    => 'gmap_vs_pmap.out';
use constant ORTHOS  => 'orthologs.tab';
use constant CI_ONLY  => 'ci_only.tab';
use constant COCCI_ONLY => 'cocci_orphans.tab';
use constant CI_ORPHANS => 'ci_orphans.tab';
#use constant LETHALS => 'lethals-current';
use constant KAKS    => 'kaks_pairwise.tab';
use constant KS      => 'kaks_pairwise.tab';
use constant KA      => 'kaks_pairwise.tab';
use constant PI      => 'popgen_10k_window.tab';
use constant SEGSITES => 'popgen_10k_window.tab';
use constant TAJIMAD => 'popgen_10k_window.tab';
use constant SNP_CI  => 'snp_ci_specific.tab';
use constant SNP_CP  => 'snp_cp_specific.tab';
use constant SNP_ALL => 'snp_all.tab';
use constant SNP_DIV => 'snp_divergences.tab';
use constant SNP_BOTH => 'snp_both_poly.tab';

#use constant WABA_BLOCKS  => 'syntenic_blocks/waba-current';
#use constant ORTHO_BLOCKS => 'syntenic_blocks/ortholog_blocks.nonoverlapping';
#use constant ORTHO_BLOCKS => 'syntenic_blocks/raw_blocks-ortholog.txt';
#use constant BLASTZ_BLOCKS => 'syntenic_blocks/raw_blocks-blastz-2k_dpmerged.txt';

use constant REPEATS  => 'repeats-rptmasker_broadrpts.tab';
#use constant PERRY   => '../repeats/perry-CeReps_CHROMOSOME_';
#use constant CLUSTER => '../families/cluster_positions_';



# TRACK COLORS - COLLECT HERE FOR EASY ACCESS
use constant BLOCKS_COLOR  => 'pink';
use constant GENES_COLOR   => 'red';
use constant ORTHOS_COLOR  => 'blue';
use constant CI_ONLY_COLOR => 'blue';
use constant CI_ORPHANS_COLOR => 'blue';
use constant COCCI_ONLY_COLOR => 'blue';
use constant GMAP_COLOR    => 'black';
use constant REPEATS_COLOR => 'green';
use constant SNP_COLOR     => 'orange';
# Font sizes...
use constant LABEL_SIZE    => '8';

# These are not being used right now
my %labels = ( 
	       genes           => 'C. immitis genes',
	       gmap            => 'gmap vs pmap',
	       orthologs       => 'orthologs / genes',
	       ci_orphans         => 'orphans / genes',
	       ci_only         => 'Ci only orphans / genes',
	       cocci_only      => 'Cocci only orphans / genes',
	       kaks            => 'average Ka/Ks',
	       ks              => 'average Ks',
	       ka              => 'average Ka',
	       mercator_blocks => 'syntenic blocks',
	       ortholog_blocks => 'Ortholog blocks',
	       repeats         => 'CocciRepeats',
	       pi              => 'pi',
	       segsites        => 'segregating sites',
	       tajimaD         => 'Tajima D (10kb windows) for C.posadasii SNPs',
	       snp_ci          => 'C.immitis Polymorphic sites',
	       snp_all         => 'All SNP sites',
	       snp_cp          => 'C.posadasii Polymorphic sites',
	       snp_div         => 'Fixed-Differences SNPs',
	       snp_both        => 'Both Polymorphic SNPs /50kb',
	       );
my %ylabels = (
	       genes           => 'genes/50kb',
	       orthologs       => 'orthologs/genes/50kb',
	       ci_orphans      => 'orphans/genes/50kb',
	       ci_only         => 'C.immitis specific/genes/50kb',
	       cocci_only      => 'Cocci specific/genes/50kb',
	       kaks            => 'mean Ka/Ks of ortholog pairs',
	       ks              => 'mean Ks of ortholog pairs',
	       ka              => 'mean Ka of ortholog pairs',
	       ortholog_blocks => 'Ortholog blocks',
	       mercator_blocks => 'Syntenic blocks',
	       repeats         => 'repetitive elements/50kb',
	       pi              => 'pi ',
	       segsites        => 'segregating sites',
	       tajimaD         => 'Tajima D (10kb windows) for C.posadasii SNPs',
	       snp_ci          => 'C.immitis Polymorphisms /50kb',
	       snp_all         => 'All SNPs /50kb',
	       snp_cp          => 'C.posadasii Polymorphisms /50kb',
	       snp_div         => 'Fixed-Differences SNPs /50kb',
	       snp_both        => 'Both Polymorphic SNPs /50kb',
	       );

my %units = (
	     genes           => 'genes/50kb',
	     gmap            => 'gmap vs pmap',
	     orthologs       => 'orthologs/genes/50kb',
	     ci_orphans      => '"orphans"/genes/50kb',
	     ci_only         => '"C.immitis specific"/genes/50kb',
	     cocci_only      => '"Cocci specific"/genes/50kb',
	     lethals         => 'lethals/genes/50kb',
	     kaks            => 'median Ka/Ks of ortholog pairs',
	     ks              => 'median Ks of ortholog pairs',
	     ka              => 'median Ka of ortholog pairs',
	     ortholog_blocks => 'Ortholog blocks',
	     mercator_blocks => 'Syntenic blocks',
	     repeats         => 'repetitive elements/50kb',
	     pi              => 'pi',
	     segsites        => 'segregating sites',
	     tajimaD         => 'Tajima D (10kb windows) for C.posadasii SNPs',
	     snp_ci          => 'C.immitis Polymorphisms /50kb',
	     snp_all         => 'All SNPs /50kb',
	     snp_cp          => 'C.posadasii Polymorphisms /50kb',
	     snp_div         => 'Fixed-Differences SNPs /50kb',
	     snp_both        => 'Both Polymorphic SNPs /50kb',
	     );

# IMAGE CONSTANTS
use constant TOP           => 35;
use constant TRACK_HEIGHT  => 80;
use constant TRACK_SPACE   => 40;
use constant TOTAL_TRACKS  => 16;

# OLD VALUES FOR TOP ALIGNED LABEL
use constant TRACK_LEFT    => 40;
use constant WIDTH         => 1024;

# Values for left aligned labels
#use constant WIDTH         => 1300;
#use constant TRACK_LEFT    => 200;

# Absolute chromosome lengths

my %CHROMS;
{
    my $dbfile = '/project/coccidiodies/variation/genomes/coccidioides_immitis_rs_2/coccidioides_immitis_rs_2.fasta';
    my $in = Bio::SeqIO->new(-format => 'fasta',-file=>$dbfile);
    my $i = 0;
    while( my $seq = $in->next_seq ) {
	my $id = $seq->display_id;
	$id =~ s/cimm(_RS_\d+)?\.(\d+)//;
	$CHROMS{$id}= [$i++, $seq->length];
    }
}

my ($largest_chrom) = sort { $b<=>$a } map {$_->[1]} values %CHROMS;  # inefficient sort

my $CONFIG = {};

# Create a bunch of objects containing the parsed data...
my $ci_genes = Windows->new('genes');
$ci_genes->parse( File::Spec->catfile($DIR,GENES),
		  qw(chrom chrom_start));
#$ci_genes->print_info('ci_genes');

my $repeats = Windows->new('repeats');
$repeats->parse(File::Spec->catfile($DIR,REPEATS),
		qw(chrom start));

my $orthos = Windows->new('orthologs');
$orthos->parse(File::Spec->catfile($DIR,ORTHOS),
	       qw(src src_start));
$orthos->normalize($ci_genes);
#$orthos->print_info('orthologs','normalized');

my $ci_orphans = Windows->new('ci_orphans'); # ci orphans
$ci_orphans->parse(File::Spec->catfile($DIR,CI_ORPHANS),
		   qw(chrom chrom_start));
$ci_orphans->normalize($ci_genes);

my $ci_only = Windows->new('ci_only'); # ci orphans
$ci_only->parse(File::Spec->catfile($DIR,CI_ONLY),
		qw(chrom chrom_start));
$ci_only->normalize($ci_genes);

my $cocci_only = Windows->new('cocci_only'); # ci orphans
$cocci_only->parse(File::Spec->catfile($DIR,COCCI_ONLY),
		   qw(chrom chrom_start));
$cocci_only->normalize($ci_genes);

my $kaks = Windows->new('kaks');
$kaks->parse(File::Spec->catfile($DIR,KAKS),
	     qw(chrom start kaks));
#$kaks->print_info('kaks','average');

my $ks = Windows->new('ks');
$ks->parse(File::Spec->catfile($DIR,KS),
	   qw(chrom start ks));

my $ka = Windows->new('ka');
$ka->parse(File::Spec->catfile($DIR,KA),
	   qw(chrom start ka));

if( 0 ) {
    my $pi = Windows->new('pi');
    $pi->parse(File::Spec->catfile($DIR,PI),
	       qw(LINKAGE_GROUP CHROM_START pi.CP));
    
    my $segsites = Windows->new('segsites');
    $segsites->parse(File::Spec->catfile($DIR,SEGSITES),
		     qw(LINKAGE_GROUP CHROM_START seg_sites.CP));
}
my $tajimaD = Windows->new('tajimaD');
$tajimaD->parse(File::Spec->catfile($DIR,TAJIMAD),
		qw(LINKAGE_GROUP CHROM_START tajima_D.CP));

my $snp_ci = Windows->new('snp_ci');
$snp_ci->parse(File::Spec->catfile($DIR,SNP_CI),
		qw(src src_start));

my $snp_cp = Windows->new('snp_cp');
$snp_cp->parse(File::Spec->catfile($DIR,SNP_CP),
		qw(src src_start));

my $snp_div = Windows->new('snp_div');
$snp_div->parse(File::Spec->catfile($DIR,SNP_DIV),
		qw(src src_start));

my $snp_all = Windows->new('snp_all');
$snp_all->parse(File::Spec->catfile($DIR,SNP_ALL),
		qw(src src_start));

my $snp_both = Windows->new('snp_both');
$snp_both->parse(File::Spec->catfile($DIR,SNP_BOTH),
		qw(src src_start));

# Just make the settings has global so I don't have to worry about
# passing it around 'cuz that just sucks
my $settings;

for my $chrom (sort { $CHROMS{$a}->[0] <=> $CHROMS{$b}->[0] } 
	       keys %CHROMS) {
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

#  plot_blocks($waba,BLOCKS_COLOR);
  plot('genes',$ci_genes,'total',GENES_COLOR);
  plot('repeats',$repeats,'total',REPEATS_COLOR);
  plot('orthologs',$orthos,'normalized',ORTHOS_COLOR);

  plot('ci_orphans',$ci_orphans,'normalized',CI_ORPHANS_COLOR);
  plot('ci_only',$ci_only,'normalized',CI_ONLY_COLOR);
  plot('cocci_only',$cocci_only,'normalized',COCCI_ONLY_COLOR);
  plot('snp_all',$snp_all,'total',SNP_COLOR);
  plot('snp_ci',$snp_ci,'total',SNP_COLOR);
  plot('snp_cp',$snp_cp,'total',SNP_COLOR);
  plot('snp_div',$snp_div,'total',SNP_COLOR);
  plot('snp_both',$snp_both,'total',SNP_COLOR);

  plot_kaks('kaks',$kaks);
  plot_kaks('ka',$ka);
  plot_kaks('ks',$ks);
  plot_posneg_numeric('tajimaD',$tajimaD);


  # Draw some header information and the xscale, which is
  # always in megabases
  my $width   = $CHROMS{$chrom}[1];
  my $scaled = TRACK_LEFT + ($xscale * $width);
  my $height = TOTAL_TRACKS * (TRACK_HEIGHT + TRACK_SPACE) + 15;
#  my $header = 'CHROMOSOME ' . $chrom;
  my $header = $chrom;
  my $footer = "megabase pairs";

  if ($GD) {
    # Place the header on the right side of the image
    $img->string(gdLargeFont,$scaled - (length($header) * gdLargeFont->width) - 2,
		 5,$header,$settings->{black});
    # ...or place it on the left edge...
    #$img->string(gdGiantFont,TRACK_LEFT,
    #	       5,$header,$settings->{black});

    $img->string(gdMediumBoldFont,$scaled/2-(length($footer)/2),
		 $height - (gdMediumBoldFont->height) - 2,$footer,$settings->{black});
    open OUT,">$chrom.png";
    print OUT $img->png;
  } else {
    $img->text(
	       style=> {
			'font' => 'Helvetica',
			'font-size'  => 40,
			'font-style' => 'bold'
		       },
	       id=>"$header",
	       x=>TRACK_LEFT,
	       y=>2)->cdata($header);

    $img->text(
	       style=> {
			'font' => 'Helvetica',
			'font-size' => 14,
		       },
	       id=>"Megabase pairs",
	       x=>$scaled/2-(length($footer)/2),
	       y=>$height)->cdata($footer);
    open OUT,">$chrom.svg";
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

  my $yscale;
  if ($ymax) {
    $yscale = TRACK_HEIGHT / $ymax;
  } else {
    $yscale = $obj->calc_y_scale($flag,TRACK_HEIGHT);
    $ymax = $obj->{max_y_total};
  }

  my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;

  # Create a grouping for this element for SVG-generated images

  my $bins = $obj->fetch_bins($chrom);
  foreach my $bin (@$bins) {
    my $total = $obj->$flag($chrom,$bin);
    my $left  = ($bin * $xscale) + TRACK_LEFT;
    my $right = $left;
    my $top   = $track_baseline + (TRACK_HEIGHT - ($total * $yscale));
    $top = ($top < $track_baseline) ? $track_baseline : $top;
    my $bottom = $track_baseline + TRACK_HEIGHT;
    if ($GD) {
      $img->rectangle($left,$top,$right,$bottom,$color);
    } else {
      $img->rectangle(x=>$left,y=>$top,
		      width  =>$right-1-$left,
		      height =>$bottom-$top,
		      id     =>"$label-$bin",
		      stroke => $color,
		      fill   => $color);
    }

    # print STDERR join("\t",$left,$top,$right,$bottom),"\n";
  }

  # Draw rectangles for each
  draw_yticks($track_baseline,$yscale,'5',$obj->{ylabel});
  draw_bounding($track_baseline,$obj);

  $CONFIG->{count}++;
  return;
}

# Draw y-scale ticks...
# This should just take a ymax, then calculate the ticks (evenly rounded) from that.
sub draw_yticks {
  my ($track_baseline,$yscale,$total_ticks,$ylabel,$sigfigs_count,
      $pos_neg) = @_;
  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count};
  my $black = $settings->{black};
  $sigfigs_count ||= 1;

  # Draw y-scale ticks...
  my $tick_left  = TRACK_LEFT;
  my $tick_right = TRACK_LEFT + 4;

  # It's printing the values in backwards order
  my $interval = TRACK_HEIGHT / $total_ticks;
  my $sigfigs = "%.".$sigfigs_count."f";  # how many significant figs for the ylabel
                    # to avoid the ylabel from not being specific
  if ( 0.20 >= ( $interval * $total_ticks) / $yscale ) {
    $sigfigs = "%.2f";
  }
  for (my $i=0; $i <= $total_ticks;$i++) {
    my $top = $track_baseline + TRACK_HEIGHT - ($i * $interval);
    my $label = ($i * $interval) / $yscale;
    # do some conditional aspect of formatting
#    my $formatted_label = sprintf(".1f",$label);
    my $formatted_label = sprintf($sigfigs,$label);
    $label ||= '0';

    if ($GD) {
      $img->line($tick_left,$top,$tick_right,$top,$settings->{black});
      $img->string(gdTinyFont,
		   $tick_left- ( length($formatted_label) * gdTinyFont->width) - 4,
		   ($top - gdTinyFont->height),$formatted_label,$settings->{black});
    } else {
      $img->line(x1 => $tick_left, y1 => $top,
		 x2 => $tick_right,y2 => $top,
		 id => "$ylabel-$formatted_label ytick",
		 stroke => $settings->{black},
		 fill   => $settings->{black});

      $img->text(
		 style=> {
			  'font' => 'Helvetica',
			  'font-size' => LABEL_SIZE,
			 },
		 id=>"$ylabel-$formatted_label",
		 x=>$tick_left-length($formatted_label) - 22,
		 y=>$top+3)->cdata($formatted_label);
    }
  }
  if( $pos_neg ) {
      for (my $i=1; $i <= $total_ticks;$i++) {
	  my $top = $track_baseline + TRACK_HEIGHT + ($i * $interval);
	  my $label = ($i * $interval) / $yscale;
	  # do some conditional aspect of formatting
	  my $formatted_label = sprintf("-".$sigfigs,$label);
	  $label ||= '0';

	  if ($GD) {
	      $img->line($tick_left,$top,$tick_right,$top,$settings->{black});
	      $img->string(gdTinyFont,
			   $tick_left- ( length($formatted_label) * gdTinyFont->width) - 4,
			   ($top - gdTinyFont->height),$formatted_label,$settings->{black});
	  } else {
	      $img->line(x1 => $tick_left, y1 => $top,
			 x2 => $tick_right,y2 => $top,
			 id => "$ylabel-$formatted_label ytick",
			 stroke => $settings->{black},
			 fill   => $settings->{black});

	      $img->text(
			 style=> {
			     'font' => 'Helvetica',
			     'font-size' => LABEL_SIZE,
			 },
			 id=>"$ylabel-$formatted_label",
			 x=>$tick_left-length($formatted_label) - 22,
			 y=>$top+3)->cdata($formatted_label);
	  }
      }
  }
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
  foreach my $gene (@sorted) {
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
    $img->string(gdTinyFont,$tick_right+4,$track_baseline,,# - ($fontheight/2),
		 '25',$settings->{black});
    # The labels still might be shifted up inappropriately
    $img->line($tick_left,
	       $track_baseline + (TRACK_HEIGHT/2),$tick_right,
	       $track_baseline + (TRACK_HEIGHT/2),
	       $settings->{black});
    $img->string(gdTinyFont,$tick_right+4,$track_baseline + (TRACK_HEIGHT/2),# - ($fontheight/2),
		 '0',$settings->{black});
    $img->line($tick_left,
	       $track_baseline + (TRACK_HEIGHT),$tick_right,
	       $track_baseline + (TRACK_HEIGHT),
	       $settings->{black});
    $img->string(gdTinyFont,$tick_right+4,$track_baseline + (TRACK_HEIGHT),# - ($fontheight/2),
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
			'font' => 'Helvetica',
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
			'font' => 'Helvetica',
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
			'font' => 'Helvetica',
			'font-size' => LABEL_SIZE,
		       },
	       id=>"gmap -25 CM",
	       x=>$tick_right + 4,
	       y=>$track_baseline + (TRACK_HEIGHT)+2)->cdata('-25');

    # Draw the unit label on the right side of the chart
    $img->text(
	       style=> {
			'font' => 'Helvetica',
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
	foreach my $point (@{$obj->{nonsliding}->{$chrom}}) {
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

    draw_yticks($track_baseline,$yscale,'5',$obj->{ylabel}, 
		1, # sig figs are 1
		1, # two-paned pos/neg plot
		);
    draw_bounding($track_baseline,$obj,1); # the third parameter indicates 
                                           # this is a pos and neg pane
    $CONFIG->{count}++;
    return;
}

sub plot_numeric { 
    my ($label,$obj) = @_;
    my $PLOT_BY = 'average';
    my $chrom  = $CONFIG->{chrom};
    my $xscale = $CONFIG->{xscale};
    my $img    = $CONFIG->{img};
    my $count  = $CONFIG->{count};
    
    my $color = $settings->{red};
    my $black = $settings->{red};
    
    # Calculate my own yscale
    my $range;
    if ($label eq 'pi') {
	$range = 0.06;
    } elsif ($label eq 'segsites') {
	$range = 100;
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
	foreach my $point (@{$obj->{nonsliding}->{$chrom}}) {
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

    draw_yticks($track_baseline,$yscale,'5',
		$obj->{ylabel}, 
		$label =~ /^pi\./ ? 2 : 1, # sig figs calculation
		0, # single pane
		);
    draw_bounding($track_baseline,$obj);
    $CONFIG->{count}++;
    return;
}

sub plot_kaks {
  my ($label,$obj) = @_;
  my $PLOT_BY = 'average';
  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count};

  my $color = $settings->{black};
  my $black = $settings->{black};
  # OLD approach using the max_y_averages...
  # That is, I calculate the values on the fly
  # my $yscale = $obj->calc_y_scale($PLOT_BY,TRACK_HEIGHT);
  # print STDERR $yscale,"\n";
  # Calculate my own yscale
  my $range;
  if ($label eq 'kaks') {
      # If we are binning, the averages are always very, very low.
      # set the range more appropriately
      #$range  = 0.5;
      #$range = 1.0;
      #$range = 1.2;
      $range = 2.0;
  } elsif ($label eq 'ks') {
      $range = 0.3;
  } elsif ( $label eq 'ka' ) {
      $range = 0.1;
  } else {
      $range = 1.0;
  }

  # The baseline for KaKs plots will be 0, centered in the middle of the plot
  # with equidistant spacing on both sides...
  my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;
  my $yscale = TRACK_HEIGHT / $range;

  # This is the old approach - plotting averages in bins across the chromosome
  my $already_plotted;
  if ( ! $already_plotted) {
    # This plots a seperate point for every ortholog pair
    foreach my $point (@{$obj->{nonsliding}->{$chrom}}) {
	my ($pos,$val) = @{$point};
	#my $fpoint = sprintf($sigfigs,$val);

	my $left  = ($pos * $xscale) + TRACK_LEFT;
	my $top   = ($track_baseline + TRACK_HEIGHT) - ($val * $yscale);
	if ( $top < $track_baseline ) {
	    $top = $track_baseline;
	}
      if ($GD) {
	  #$img->setPixel($left,$top,$color);
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
  # for KA/KS
  if ($label eq 'kaks') {
    my $width  = $CHROMS{$chrom}[1];
    my $right  = TRACK_LEFT + ($xscale * $width) + 2;
    my $location = ($track_baseline + TRACK_HEIGHT) - (1 * $yscale);
    if ($label eq 'kaks') {
      if ($GD) {
	$img->line(TRACK_LEFT,$location,$right,$location,$settings->{aqua});
      } else {
	$img->line(x1=>TRACK_LEFT,
		   y1=>$location,
		   x2=>$right,
		   y2=>$location,
		   id=>"$label-significance divider",
		   stroke => $settings->{aqua},
		   fill   => $settings->{aqua});
      }
    }
  }

  draw_yticks($track_baseline,
	      $yscale,'5',$obj->{ylabel},$range < 0.5 ? "2" : '1');
  draw_bounding($track_baseline,$obj);
  $CONFIG->{count}++;
  return;
}


sub plot_blocks {
  my ($obj,$icolor) = @_;
  my $color = $settings->{$icolor};
  my $black = $settings->{black};

  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count};

  my $track_baseline = TOP + ((TRACK_HEIGHT + TRACK_SPACE) * $count) + TRACK_SPACE / 2;

  my @blocks = @{$obj->{$chrom}};
  my $total;
  foreach my $block (@blocks) {
    $total++;
    my ($start,$stop) = @$block;
    ($start,$stop) = ($stop,$start) if ($start > $stop);
    my $left   = ($start * $xscale) + TRACK_LEFT;
    my $right  = ($stop * $xscale) + TRACK_LEFT;
    my $top    = $track_baseline;
    my $bottom = $track_baseline + TRACK_HEIGHT;
    if ($GD) {
      $img->filledRectangle($left,$top+1,$right-1,$bottom-1,$color);
    } else {
      $img->rectangle(x=>$left,y=>$top,
		      width  =>$right-1-$left,
		      height =>$bottom-$top,
		      id     =>"waba_blocks-$total-" . $start .'-' . $stop,
		      stroke =>$color,
		      fill   =>$color),
    }
    # print STDERR join("\t",$left-1,$top,$right+1,$bottom),"\n";
  }

  draw_bounding($track_baseline,$obj);
  $CONFIG->{count}++;
  return;
}

# Draw a bounding box for this track.
sub draw_bounding {
  my ($top,$obj,$two_pane) = @_;
  my $track_label = $obj->{ylabel};
  my $units       = $obj->{units};
  my $chrom  = $CONFIG->{chrom};
  my $xscale = $CONFIG->{xscale};
  my $img    = $CONFIG->{img};
  my $count  = $CONFIG->{count};

  my $width  = $CHROMS{$chrom}[1];
  my $left   = TRACK_LEFT;
  my $right  = TRACK_LEFT + ($xscale * $width);
  my $bottom = $top + TRACK_HEIGHT;
  if( $two_pane ) {
      $bottom += TRACK_HEIGHT;
  }

  if ($GD) {
    $img->rectangle($left,$top,$right,$bottom,$settings->{black});
  } else {
    $img->rectangle(x     =>$left,
		    y     =>$top,
		    width =>$right-$left,
		    height=>$bottom-$top,
		    id=>$track_label,
		    stroke=>$settings->{black},
		    'fill-opacity'=>0);
  }

  # Draw xscale tickmarks every million basepairs
  my $tick_top     = $top + TRACK_HEIGHT - 4;
  my $tick_bottom  = $top + TRACK_HEIGHT;

  if( $two_pane ) { 
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
      $img->string(gdTinyFont,
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
			  'font' => 'Helvetica',
			  'font-size' => LABEL_SIZE,
			  'color'=> $settings->{black},
			 },
		 id=>"$track_label-$label MBp",
		 x=>$left-3,
		 y=>$tick_bottom + 10)->cdata($label);
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
			'font' => 'Helvetica',
			'font-size' => 14,
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
		  purple => $img->colorAllocate(200,100,200),
		  grey   => $img->colorAllocate(230,230,230),
		  pink   => $img->colorAllocate(204,000,204),
		  sea    => $img->colorAllocate(000,102,102),
		  orange => $img->colorAllocate(255,153,000),
		  # Good colors
		  red    => $img->colorAllocate(255,0,0),
		  blue   => $img->colorAllocate(000,000,255),
		  green  => $img->colorAllocate(000,255,000),
		  aqua   => $img->colorAllocate(000,255,204),
		  
		  fontwidth  => gdMediumBoldFont->width,
		  fontheight => gdMediumBoldFont->height,
		};
  } else {
    $settings = { white  => 'rgb(255,255,255)',
		  black  => 'rgb(0,0,0)',
		  purple => 'rgb(200,100,200)',
		  grey   => 'rgb(230,230,230)',
		  pink   => 'rgb(204,000,204)',
		  sea    => 'rgb(000,102,102)',
		  orange => 'rgb(255,153,000)',
		  # Good colors
		  red    => 'rgb(255,0,0)',
		  blue   => 'rgb(000,000,255)',
		  green  => 'rgb(000,255,000)',
		  aqua   => 'rgb(000,255,204)',
		  #		  fontwidth  => gdMediumBoldFont->width,
		  #		  fontheight => gdMediumBoldFont->height,
		};
  }
  return $settings;
}


package Windows;

use List::Util qw(sum max);
# constants for sliding windows
use constant WINDOW        => 50_000;
use constant STEP          => 10_000;

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
  my ($self,$file,$group_by,$bin_by,$save_by) = @_;
  my $cols = fetch_columns($file);
  for my $col ( $group_by, $bin_by, $save_by) {
      if( defined $col && ! exists $cols->{$bin_by}) {
	  die("cannot find column $bin_by in $file\n");
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
      print STDERR "Positional value out of range for chromosome...$bin_val\t",
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

	# HORRENDOUS KLUDGE!
	# ignore ka/ks values of 0.  sumthin's wrong.
	# Need to clip outlier values (ks > 5 and Ka > 3);
	if ($save_by eq 'kaks' || $save_by eq 'ks' || $save_by eq 'ka') {
	    my $kaks = $fields[$cols->{kaks}];
	    my $ks   = $fields[$cols->{ks}];
	    my $ka   = $fields[$cols->{ka}];
	    next if ($ks > 1.5 || $ka > 4 );
	    #if ($value == 0) {
	    #	print STDERR "WARNING: Ka/Ks value out of range : $value\n";
	    #	print STDERR "\tPositional information\t",join("\t",$bin_val,$value,$group_val),"\n";
	    #	next;
	    #}
	    # Save every point - I might just want to plot the raw data itself
	    push (@{$self->{nonsliding}->{$group_val}},[$bin_val,$value]);
	} elsif( $save_by =~ /^(pi|seg_sites|tajima_D)/ ) {
	    push (@{$self->{nonsliding}->{$group_val}},[$bin_val,$value]);
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
  my ($self,$file) = @_;
  print STDERR "parsing: $file ...\n";
  
  my $positions = {};
  open IN,$file or die "$! $file\n";
  while (<IN>) {
    chomp;
    # Skip comments
    next if (/^\#/);
    my($chrom,$start,$end,$strand,$target) = split "\t";
    # Fetch the position of the bin_by and group_by
    # columns in the fields array
    push (@{$self->{$chrom}},[$start,$end]);
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
    foreach my $group_value (keys %$positions) {
	# Get the maximal limit (ie the highest scoring feature)
	my @positions = sort { $a->[0] <=> $b->[0] } @{$positions->{$group_value}};
	foreach my $temp (@positions) {
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
#    print STDERR "\t",$i,"\t",$bin,"\t",$value,,"\t",$fstart,"\n";
    }
    return;
}


# Calculate an average for each bin or a total value.
sub calc_averages {
    my $self = shift;
    warn("   calculating average...\n");
    my ($max_average,$max_total);
    for my $group_by ($self->fetch_groups()) {	
	my $bins = $self->fetch_bins($group_by);
	for my $bin (@$bins) {
	    my $all = $self->fetch_values($group_by,$bin);
	    my $sum = &sum (@$all);
	    my $avg = $sum / scalar @$all;
# print STDERR $bin,"\t",$group_by,"\t",$sum,"\t",scalar @$all,"\t",$avg,"\n";
#	    warn( $self->{label}, " AVERAGE > 1000\n",join("\t\n",@$all),"\n") if $avg > 1000;
	    $self->{ max_y_average} = (! defined $self->{max_y_average} || 
				       $avg > $self->{max_y_average}) ? 
		$avg : $self->{max_y_average};
	    $max_total   = (! defined $max_total || 
			    scalar @$all > $max_total) ? scalar @$all : 
			    $max_total;

	    $self->{groups}->{$group_by}->{$bin}->{average} = $avg;
	    $self->{groups}->{$group_by}->{$bin}->{total} = scalar @$all;

	    $self->{max_y_total}   = $max_total;
	}
    }
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
    foreach (@cols) {
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
    foreach my $group ($numerator->fetch_groups()) {
	my $bins = $numerator->fetch_bins($group);
	foreach my $bin (@$bins) {
	    my $total = $numerator->total($group,$bin);
	    my $denom_total = $denominator->total($group,$bin);
	    my $normalized = eval { $total / $denom_total };
	    $normalized ||= '0';

	    $max = (! defined $max || 
		    $normalized > $max) ? $normalized : $max;
	    $numerator->{groups}->{$group}->{$bin}->{normalized} = $normalized;
	    $numerator->{max_y_normalized} = $max;
	}
    }
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
  foreach my $group ($self->fetch_groups()) {

    # Print out some header information
    print OUT "#chromosome=$group\n";
    print OUT "#field_content=$field\n";
    print OUT "#max_y_average=" . $self->{max_y_average} . "\n";
    print OUT "#max_y_total=" . $self->{max_y_total} . "\n";
    print OUT "#max_y_normalized=" . $self->{max_y_normalized} . "\n";
    my $bins = $self->fetch_bins($group);
    foreach my $bin (@$bins) {
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
  my $scale = $height / $max;
  return $scale;
}



1;
