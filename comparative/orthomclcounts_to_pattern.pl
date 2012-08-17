#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $style_float = 
'    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.min.js"></script>

    <script type="text/javascript">
        function UpdateTableHeaders() {
            $("div.divTableWithFloatingHeader").each(function() {
                var originalHeaderRow = $(".tableFloatingHeaderOriginal", this);
                var floatingHeaderRow = $(".tableFloatingHeader", this);
                var offset = $(this).offset();
                var scrollTop = $(window).scrollTop();
                if ((scrollTop > offset.top) && (scrollTop < offset.top + $(this).height())) {
                    floatingHeaderRow.css("visibility", "visible");
                    floatingHeaderRow.css("top", Math.min(scrollTop - offset.top, $(this).height() - floatingHeaderRow.height()) + "px");

                    // Copy cell widths from original header
                    $("th", floatingHeaderRow).each(function(index) {
                        var cellWidth = $("th", originalHeaderRow).eq(index).css(\'width\');
                        $(this).css(\'width\', cellWidth);
                    });

                    // Copy row width from whole table
                    floatingHeaderRow.css("width", $(this).css("width"));
                }
                else {
                    floatingHeaderRow.css("visibility", "hidden");
                    floatingHeaderRow.css("top", "0px");
                }
            });
        }

        $(document).ready(function() {
            $("table.tableWithFloatingHeader").each(function() {
                $(this).wrap("<div class=\"divTableWithFloatingHeader\" style=\"position:relative\"></div>");

                var originalHeaderRow = $("tr:first", this)
                originalHeaderRow.before(originalHeaderRow.clone());
                var clonedHeaderRow = $("tr:first", this)

                clonedHeaderRow.addClass("tableFloatingHeader");
                clonedHeaderRow.css("position", "absolute");
                clonedHeaderRow.css("top", "0px");
                clonedHeaderRow.css("left", $(this).css("margin-left"));
                clonedHeaderRow.css("visibility", "hidden");

                originalHeaderRow.addClass("tableFloatingHeaderOriginal");
            });
            UpdateTableHeaders();
            $(window).scroll(UpdateTableHeaders);
            $(window).resize(UpdateTableHeaders);
        });
    </script>

<style>
<!--
th {
    background-color: lightgrey;
    border: 1px solid black;
}
td {
    border: 1px solid black;
}
-->
</style>
';

my $hmmurl = "http://pfam.sanger.ac.uk/family/%s";
my $seqfolder = "../../%s.orthomcl_seqs";

GetOptions("hmm|url:s" => \$hmmurl,
	   'folder|seq:s' => \$seqfolder);
my $file = shift;
my $base = $file;
if( $file =~ /(\S+)\.(report|dat|tab|txt)$/ ) {
	$base = $1;	
}
my ($idbase) = split(/\./,$base); # some assumtions about the naming procedure
$seqfolder = sprintf($seqfolder,$idbase); 

mkdir("$base.family_by_pattern.d");
open(my $fh => $file) || die "$file: $!";
my $header = <$fh>;
my ($name,@strains) = split(/\t/,$header);
pop @strains;
my %patterns;
my %fams;
my %pat2name;
my $counter = 0;
while (<$fh>) {
  my @row= split(/\t/,$_);
  my ($family,@counts) = @row;
  my $desc = pop @counts;
  # add URL to desc
  my @newdesc;
  for my $domain ( split(/\s*,\s*/,$desc ) ) {
   if( $domain =~ /(\S+)\s+\((\d+)\)/ ) {
     my $url = sprintf($hmmurl,$1);
     $domain = sprintf("<a target=pfam href=\"%s\">%s</a> (%d)",$url,$1,$2);
   }
   push @newdesc, $domain;
  }
  $desc = join(",", @newdesc);
  my $i = 0;
  my (@pattern,@name_pattern);
  for my $c ( @counts ) {
   if( $c ) {
      push @name_pattern, $strains[$i];
    }
    push @pattern, ( $c > 0 )  ? 1 : 0;
    $i++;
  }
  if( ! exists $pat2name{join(",",@pattern)} ) { 
     $pat2name{join(",",@pattern)} = $counter++;
  }
#join(",", @name_pattern);
  $patterns{join(",",@pattern)}++;
  push @{$fams{join(",",@pattern)}}, [sprintf("<a target=seqfam href=\"%s/%s.pep.fa\">%s</a>",$seqfolder,$family,$family), @counts, $desc];
}

open(my $ofh => ">$base.patterns.tab") || die $!;
print $ofh join("\t", qw(COUNT),@strains),"\n";
for my $p ( sort { $patterns{$b} <=> $patterns{$a} } keys %patterns ) {
  print $ofh join("\t", $patterns{$p}, split(/,/,$p)),"\n";
}

open($ofh => ">$base.patterns.html") || die $!;
print $ofh "<HTML><HEAD><title>$base family report</title>\n$style_float</HEAD>\n<BODY BGCOLOR=\"WHITE\">\n";
print $ofh "<table class=\"tableWithFloatingHeader\" style=\"border: 3px\">\n<thead><tr>\n  ";
for my $nm ( qw(COUNT), @strains) {
 printf $ofh "<th>%s</th>",$nm;
}
print $ofh "</tr></thead>\n<tbody>\n";
for my $p ( sort { $patterns{$b} <=> $patterns{$a} } keys %patterns ) {
  my $name = $pat2name{$p};
  print $ofh "<tr>\n  "; 
  printf $ofh "<td><a href=\"%s\">%s</td>\n","$base.family_by_pattern.d/$name.report.html",$patterns{$p};
  for my $cell ( split(/,/,$p) ) {
      printf $ofh " <td BGCOLOR=\"%s\" />\n",$cell ? 'BLACK' : 'WHITE';
  }
  print $ofh "</tr>\n";
}
print $ofh "</tbody>\n</table>\n</BODY><FOOT></FOOT></HTML>\n";

for my $pat ( keys %fams ) {
 my $name = $pat2name{$pat};
 open($ofh => ">$base.family_by_pattern.d/$name.report.html") || die "$!";
 print $ofh "<HTML><HEAD><title>families with pattern $name</title>\n$style_float</HEAD><BODY BGCOLOR=\"WHITE\">\n";
 print $ofh "<table class=\"tableWithFloatingHeader\" style=\"border: 3px\">\n<thead><tr>\n  ";
 for my $nm ( qw(FAMILY), @strains,'PFAM Domains') {
  printf $ofh "<th>%s</th>",$nm;
 }
 print $ofh "</tr></thead>\n<tbody>\n";
 for my $d ( @{$fams{$pat}} ) {
  print $ofh "<tr>\n ";
  for my $col ( @$d ) {
   printf $ofh "  <td>%s</td>\n",$col;
  }
  print $ofh "</tr>\n";
 }
 print $ofh "</tbody>\n</table>\n</BODY><FOOT></FOOT></HTML>\n";
}
