#!/usr/bin/perl -w
# This script read in original network file and convert to integer ID
# usage: perl src/parse_net.pl <original network file>
# output: 
# # 1. <original network file>_id : e.g: test_data/Yeast_TF_net.txt_id
# # 2. <original network file>_id_net : e.g: test_data/Yeast_TF_net.txt_id_net

if($#ARGV < 0){
	print "usage: perl src/parse_net.pl <original network file>\n";
	die;
}
$input = $ARGV[0]; ## input file
$out1  = $input."_id";
$out2  = $input."_id_net";

open I,$input or die $!;
open O1,">$out1" or die $!;
open O2,">$out2" or die $!;
$k = 1;$del=0;
while(<I>){
	chomp;
	($a,$b,$w) = split " ";
	($a,$b)=sort($a,$b);
	$mark=$a."_".$b;
	if($e{$mark}){
		$del++;
		next;
	}
	$e{$mark}=1;
	unless($id{$a}){
		$id{$a} = $k;
		print O1 $a."\t".$k."\n";	
		$k++;
	}
	unless($id{$b}){
		$id{$b} = $k;
		print O1 $b."\t".$k."\n";	
		$k++;
	}
	$aa = $id{$a};
	$bb = $id{$b};
	print O2 $aa."\t".$bb."\t".$w."\n";
}
print "Delete duplicate edges: $del\n";
close I;
close O1;
close O2;
