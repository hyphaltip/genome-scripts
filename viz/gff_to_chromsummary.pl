#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqFeature::Generic;
use File::Spec;
use Bio::Graphics;

my $WIDTH = 1000;

my $dir = 'graphics';
my $type_filter;
GetOptions('d|dir:s' => \$dir,
	   't|type:s' => \$type_filter,
	   );
my $db = shift;
my $gff = shift;

mkdir($dir) unless -d $dir;
my $dbh = Bio::DB::Fasta->new($db);
open(my $fh => $gff) || die $!;
my %dat;
while(<$fh>) {
    next if /^\#/;
    my ($chrom,$src,$type,$start,$stop,$score,$strand,$frame,
	$group) = split;
    next if( defined $type_filter && $type_filter ne $type);

    my %group = map { split(/=/,$_) } split(/;/,$group);
    push @{$dat{$chrom}}, Bio::SeqFeature::Generic->new(-start => $start,
							-end   => $stop,
							-primary_tag => $type,
							-source_tag  => $src,
							-tag => { Note => $group{Note} },
							-display_name => $group{Name});
}
my $limit = $type_filter || 'all';
for my $chrom ( keys %dat ) {    
    my $ofname = File::Spec->catfile($dir,"$chrom.$limit.png");
    my $chrom_length = $dbh->length($chrom);
    my $panel = Bio::Graphics::Panel->new(-length=> $chrom_length,
					  -width  => $WIDTH,
					  -pad_left => 0,
					  -pad_right => 0);
    
    my $top_length = Bio::SeqFeature::Generic->new
	(-start => 1,
	 -end   => $chrom_length,
	 -display_name => $chrom);
    
    $panel->add_track($top_length,
		      -glyph   => 'arrow',
		      -tick    => 2,
		      -fgcolor => 'black',
		      -double  => 1,
		      -label   => 1,
		      );
    $panel->add_track($dat{$chrom},
		      -glyph => 'transcript',
		      -label => 0,
		      -bump  => 0,
		      -bgcolor => 'blue',
		      -fgcolor => 'blue',
		      -height => 5,
#		      -description => \&desc,
		      -font2color => 'red');    
    open(my $ofh => ">$ofname") || die $!;
    binmode($fh);
    print $ofh $panel->png;
}

sub desc {
    my $ft = shift;
    my ($n) = $ft->get_tag_values('Note');
    $n;
}
