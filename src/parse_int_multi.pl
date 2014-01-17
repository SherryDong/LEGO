#!/usr/bin/perl -w
# This script read in interesting gene list file and convert to integer ID

if($#ARGV < 1){
	print "usage: perl src/parse_int_multi.pl <interesting gene list file> <gene Id file> \n";
	die;
}
$input = $ARGV[0]; ## input gene set file
$id_f  = $ARGV[1]; ## input gene Id file
$out1  = $input."_id";
$out2  = $input."_id_2";
##
open ID,$id_f or die $!;
while(<ID>){
	chomp;
	($g,$i) = split "\t";
	$id{$g} = $i;
}
close ID;
##
open I,$input or die $!;
open O1,">$out1" or die $!;
open O2,">$out2" or die $!;
$k = 0;
while(<I>){
	chomp;
	($a,$b) = split "\t",$_;
	unless($idd{$b}){
		$k++;
		print O2 $b."\t".$k."\n";
		$idd{$b} = $k;
	}
	if($id{$a}){
		$aa = $id{$a};
		$bb = $idd{$b};
		$final{$aa."\t".$bb}=1;
	}else{
		$no{$a}=1;
	}
}
close I;
@n = sort(keys %final);
$n = join "\n",@n;
print O1 $n."\n";
close O1;
@no = keys %no;$no=$#no+1;
if($no>0){
	print "There are $no genes from interesting gene list not found in the $id_f!\n";
}
