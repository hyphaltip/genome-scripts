#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use List::Util qw(sum);

use constant { 
    style_float => '    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.min.js"></script>

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
',
};

my $hmmurl = "http://pfam.sanger.ac.uk/family/%s";
my $ext = 'domtbl.tab';
GetOptions("hmm|url:s" => \$hmmurl);
my $dir = shift || die "must profile a dir with hmmscan output!";
opendir(DIR,$dir) || die $!;
my %domains;
my %domain_desc;
my %spall;
for my $file ( readdir(DIR) ) {
    next unless ($file =~ /\.\Q$ext\E$/);
    open(my $fh => "$dir/$file") || die "$dir/$file: $!";
    while(<$fh>) {
	next if /^\#/;
	chomp;
	my ($pfam_name, $pfam_acc, $tlen, $qname, $qacc, $qlen, $full_evalue,@rest) = split(/\s+/,$_,23);
	my ($sp) = split(/\|/,$qname);
	$domains{$pfam_name}->{$sp}++;
	$domain_desc{$pfam_name} = pop @rest unless exists $domain_desc{$pfam_name};
	$spall{$sp}++
    }
}
      
print "<HTML><HEAD><title>Pfam family report</title>\n",style_float,"</HEAD>\n<BODY BGCOLOR=\"WHITE\">\n";
print "<table class=\"tableWithFloatingHeader\" style=\"border: 3px\">\n<thead><tr>\n  ";
my @sp = sort keys %spall;
for my $nm ( qw(Pfam Desc), @sp) {
 printf "<th>%s</th>",$nm;
}
print "</tr></thead>\n<tbody>\n";

for my $domain (map { $_->[1] } 
		sort { $b->[0] <=> $a->[0] }
		map { [ (sum values %{$domains{$_}}), $_ ] }
		keys %domains ) {
  print "<tr>\n ";
  printf "<td><a target=pfam href=\"$hmmurl\">%s</a></td><td>%s</td>",$domain,$domain,$domain_desc{$domain};
  for my $col ( @sp ) {
   printf "  <td>%s</td>\n",exists $domains{$domain}->{$col} ? $domains{$domain}->{$col} : 0;
  }
  print "</tr>\n";
}
print "</tbody>\n</table>\n</BODY><FOOT></FOOT></HTML>\n";

