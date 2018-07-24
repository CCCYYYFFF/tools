#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path);

my ($in,$out);
GetOptions(
   "in:s"  =>\$in,
   "out:s"  =>\$out,
) or &USAGE;

&USAGE unless ($in);
$out||="huizong.txt";
$in = abs_path($in);
my @name = `ls $in`;
chomp(@name);
open(O,">$out");

print O "name\t\n";
foreach my $i (@name){
	my($a,$b);
	my $file = "$in/$i/Report";  
	open(I,"$file/report_QC.txt");
	print O "$i |\t";
	while(<I>){
		chomp;
		my @ord = split(/\t+/,$_);
		if($_=~/good quality/ && $_=~/no muts/){
			print O "$ord[2]";
		}
		elsif($_=~/good quality/ && $_=~/have muts/){
			my $var = `ls $file/$ord[0]|grep chr`;
			$var=~/.+\.(.+)\.svg/;
			print O "$1";
		}
		else{
			print O "$ord[1],$ord[2]";
		}
		print O "|\t";
	}
	close(I);
	print O " &&\n";
}
close(O);

sub USAGE {
 my $usage = <<"USAGE";

USAGE:
 -in input file dir
 -out out file dir
USAGE
 print $usage;
 exit;
}
