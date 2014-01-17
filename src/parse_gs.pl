#!/usr/bin/perl -w
# This script read in original gene set file and convert to integer ID
# usage: perl src/parse_gs.pl <original gene set file> <gene Id file>
# output: 
# <original gene set file>_id : e.g : test_data/Gene_Set_list.txt_id
# <original gene set file>_id_gs : e.g: test_data/Gene_Set_list.txt_id_gs

if($#ARGV < 1){
	print "usage: perl src/parse_gs.pl <original gene set file> <gene Id file> \n";
	die;
}
$input = $ARGV[0]; ## input gene set file
$id_f  = $ARGV[1]; ## input gene Id file
$out0  = $input."_gene_id";
$out1  = $input."_id";
$out2  = $input."_id_gs";
##
open ID,$id_f or die $!;
while(<ID>){
	chomp;
	($g,$i) = split "\t";
	$id{$g} = $i;
}
close ID;
$max_id = $i;
##
open I,$input or die $!;
system "cp $id_f $out0";
open O0,">>$out0" or die $!;
open O1,">$out1" or die $!;
open O2,">$out2" or die $!;
$k = 1;
while(<I>){
	chomp;
	($a,$b) = split "\t";
	unless($id{$a}){
		$no{$a} = 1;		
		$max_id ++;
		print O0 $a."\t".$max_id."\n";
		$id{$a} = $max_id;
	}
	if($id{$a}){
		unless($id{$b}){
			$id{$b} = $k;
			print O1 $b."\t".$k."\n";	
			$k++;
		}
		$aa = $id{$a};
		$bb = $id{$b};
		print O2 $aa."\t".$bb."\n";
	}
}
close I;
close O0;
close O1;
close O2;
@no = keys %no;$no=$#no+1;
if($no>0){
	print "There are $no genes not found in the network;Add them in $out0;\n";
}
