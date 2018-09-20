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
$tag||="p";
$note||="T";
$in = abs_path($in);
my @name = `find $in -name "*.svg"`;
chomp(@name);

foreach my $x (@name){
	open(IN,"$x");
	open(OUT,">$x.bak");
	my ($old,$pos);
	while(<IN>){
		if($_ =~/circle.*cx="(\d+)"/){
			$old = $_;
			$pos = $1;
			#print OUT "<ellipse cx=\"$1\" cy=\"360\" rx=\"40\" ry=\"400\" style=\"fill:blue;fill-opacity:0;stroke-opacity:1;stroke:red;stroke-width:4\" />\n";
			next;
		}elsif($_ =~/chr\d+,\d+,(\w+),(\w+),/){
			if(length($1)==1 && length($2)==1 && $tag eq "e"){
				print OUT "<ellipse cx=\"$pos\" cy=\"360\" rx=\"40\" ry=\"400\" style=\"fill:blue;fill-opacity:0;stroke-opacity:1;stroke:red;stroke-width:4\" />\n";
			}else{
				print OUT "$old";
			}
			if($note eq "F"){
				next;
			}
		}
		print OUT $_;
	}
	close(IN);
	close(OUT);
	`mv $x $x.old`;
	`mv $x.bak $x`;
}

sub USAGE {
 my $usage = <<"USAGE";

USAGE:
 -in input file dir
 -tag point(p) or ellipse(e)
 -note have note(T) or not(F)
USAGE
 print $usage;
 exit;
}
