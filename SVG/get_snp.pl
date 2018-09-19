#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);

my ($in,$tag,$note);
GetOptions(
   "in:s"  =>\$in,
   "tag:s"  =>\$tag,
   "note:s" =>\$note,
) or &USAGE;

&USAGE unless ($in);
$tag||="T";
$note||="p";
$in = abs_path($in);
open(IN,"$in");
open(OUT,">$in.bak");

while(<IN>){
	if($_ =~/circle.*cx="(\d+)"/ && $tag eq "e"){
		print OUT "<ellipse cx=\"$1\" cy=\"360\" rx=\"40\" ry=\"400\" style=\"fill:blue;fill-opacity:0;stroke-opacity:1;stroke:red;stroke-width:4\" />\n";
		next;
	}elsif($_ =~/chr\d+,/ && $note eq "F"){
		next;
	}
	print OUT $_;
}

close(IN);
close(OUT);

`mv $in $in.old`;
`mv $in.bak $in`;

sub USAGE {
 my $usage = <<"USAGE";

USAGE:
 -in input SVG file 
 -tag point(p) or ellipse(e)
 -note have note(T) or not(F)
USAGE
 print $usage;
 exit;
}
