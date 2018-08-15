#!/usr/bin/perl

#*****************************************************************************
# FileName: get_mut.pl
# Creator: Chenyongfeng <chenyongfeng@celloud.cn>
# Create Time: 2018-8-15
# Description: get the codon mut.
# CopyRight: Copyright (c) CelLoud, All rights reserved.
# Revision: V1.0
#*****************************************************************************

use strict;
use warnings;
use utf8;
use Cwd 'abs_path';
use FindBin qw($Bin);
use File::Basename qw(basename dirname);

if (@ARGV!=2) {
	die"usage : *.pl <log.txt> <Mut.txt>\n";
	exit;
}

my %codon = ("TTT"=>"F","TCT"=>"S","TAT"=>"Y","TGT"=>"C",      
	  "TTC"=>"F","TCC"=>"S","TAC"=>"Y","TGC"=>"C",      
	  "TTA"=>"L","TCA"=>"S","TAA"=>"*","TGA"=>"*",
	  "TTG"=>"L","TCG"=>"S","TAG"=>"*","TGG"=>"W",   
	  "CTT"=>"L","CCT"=>"P","CAT"=>"H","CGT"=>"R",      
	  "CTC"=>"L","CCC"=>"P","CAC"=>"H","CGC"=>"R",      
	  "CTA"=>"L","CCA"=>"P","CAA"=>"Q","CGA"=>"R",      
	  "CTG"=>"L","CCG"=>"P","CAG"=>"Q","CGG"=>"R",      
	  "ATT"=>"I","ACT"=>"T","AAT"=>"N","AGT"=>"S",      
	  "ATC"=>"I","ACC"=>"T","AAC"=>"N","AGC"=>"S",      
	  "ATA"=>"I","ACA"=>"T","AAA"=>"K","AGA"=>"R",      
	  "ATG"=>"M","ACG"=>"T","AAG"=>"K","AGG"=>"R",      
	  "GTT"=>"V","GCT"=>"A","GAT"=>"D","GGT"=>"G",      
	  "GTC"=>"V","GCC"=>"A","GAC"=>"D","GGC"=>"G",      
	  "GTA"=>"V","GCA"=>"A","GAA"=>"E","GGA"=>"G",      
	  "GTG"=>"V","GCG"=>"A","GAG"=>"E","GGG"=>"G" );

my %original_mut=(	"169"=>"I169T", 
		"173"=>"V173L",
		"180"=>"L180M",
		"181"=>"A181V|T",
		"184"=>"T184G",
		"202"=>"S202G|I",
		"204"=>"M204V|I",
		"236"=>"N236T",
		"250"=>"M250V",
		"214"=>"V214A",
		"215"=>"Q215S",
		"229"=>"L229W|V",
		"233"=>"I233V",	);

open IN,"$ARGV[0]";
my %original;
my %hash1=%original_mut;
my %codon_new;
my %codon_new2;
my %lqc;
my %jianji_old;
my %score;
my %SNP;
my %postion;
while (<IN>){
	chomp;
	if(/^Phred_seq/){next;}
	if(/^$/){next;}
	if(/nucleotides was detected/){next;}
	if(/^\#Type:/){next;}
	if (/^\d+:(\S+):\S+/){
		my @temp = split(/:/,);
		$original{$temp[0]} = $temp[2];
		$hash1{$temp[0]} = 0;
	}
	my @temp = split/\t/;
	if(exists $temp[5]){
		my @pos_three = split(/\:/,$temp[2]);
		my @jianji_three = split(//,$temp[3]);
		if ($temp[5] =~/(\w)(\d+)/){
			my $pep_pos = $2;
			$postion{$pep_pos} = $temp[2];
			$jianji_old{$pep_pos} = $temp[3];
			if(exists $temp[6]){
				$lqc{$pep_pos} = $temp[6];
			}else{
				$lqc{$pep_pos} = "good qc";
			}
			if(!exists $codon_new{$pep_pos}){
				$original{$pep_pos} = $1;
				$hash1{$pep_pos} = 0;
				$codon_new{$pep_pos}{$pos_three[0]} = $jianji_three[0];
				$codon_new{$pep_pos}{$pos_three[1]} = $jianji_three[1];
				$codon_new{$pep_pos}{$pos_three[2]} = $jianji_three[2];
			}
			if($temp[4] =~/(.) - (.)\|(.)\(.*,(.*)\)/){
				$SNP{$temp[1]} = 2;
				$score{$pep_pos} .= $temp[4].";";
				$codon_new{$pep_pos}{$temp[1]} = $2;
				$codon_new2{$pep_pos}{$temp[1]} = $3;
			}elsif($temp[4] =~/(.) - (.)/){
				if ($1 ne $2){
					$SNP{$temp[1]} = 1;
					$score{$pep_pos} .= $temp[4].";";
					$codon_new{$pep_pos}{$temp[1]} = $2;
				}
			}
		}
	}
}
close IN;
my %mutation;
open OUT,">$ARGV[1]";
foreach my $key1 (sort {$a<=>$b} keys %hash1){
	if(!exists $original{$key1}){
		print OUT "$original_mut{$key1}\tNO Find the Loci\n";
		next;
	}
	if(!exists $codon_new{$key1}){
		print OUT "$original_mut{$key1}\tNO Mut\n";
		next;
	}
	my $codon_lishi;
	my $codon_lishi2;
	my $new_wz="";
	$mutation{$key1} = "$original{$key1}$key1";
	foreach my $key2 (sort{$a<=>$b} keys %{$codon_new{$key1}}){
		$codon_lishi .= $codon_new{$key1}{$key2};
		if(exists $SNP{$key2}){
			$new_wz .= $key2.";";
			if($SNP{$key2}==2){
				$codon_lishi2 .= $codon_new2{$key1}{$key2};
			}else{
				$codon_lishi2 .= $codon_new{$key1}{$key2};
			}
		}else{
			$codon_lishi2 .= $codon_new{$key1}{$key2};
		}
	}
	$mutation{$key1} .= "$codon{$codon_lishi}";
	if($codon{$codon_lishi} ne $codon{$codon_lishi2}){
		$mutation{$key1} .= "|$codon{$codon_lishi2}";
		$codon_lishi2 =";".$codon_lishi2;
	}else{
		$codon_lishi2 = "";
	}
	print OUT "$mutation{$key1}\t$lqc{$key1}\t$key1\t$new_wz\t$postion{$key1}\t$jianji_old{$key1}-$codon_lishi$codon_lishi2\t$score{$key1}\n";
}
close OUT;
