#!/usr/bin/perl -w
# This script read in interesting gene list file and convert to integer ID

if($#ARGV < 1){
	print "usage: perl src/parse_int_multi_bg.pl <interesting gene list file> <background file> <gene Id file> <use bg filter or not>>\n";
	die;
}
$input = $ARGV[0]; ## input gene set file
$bg    = $ARGV[1]; ## input gene set file
$id_f  = $ARGV[2]; ## input gene Id file
$filter = $ARGV[3];
$out1  = $input."_id";
$out2  = $bg."_id";
$out11  = $input."_id_2";
##
open ID,$id_f or die $!;
while(<ID>){
	chomp;
	($g,$i) = split "\t";
	$id{$g} = $i;
}
close ID;
##
open I,$bg or die $!;
open O2,">$out2" or die $!;
while(<I>){
	chomp;
	@a = split "\t";
	$a = $a[0];
	if($id{$a}){
		$aa = $id{$a};
		$bg{$aa} = 1;
	}else{
		$no{$a}=1;
	}
}
close I;
@n = sort{$a<=>$b}(keys %bg);
$n = join "\n",@n;
print O2 $n."\n";
close O2;
##
open I,$input or die $!;
open O1,">$out1" or die $!;
open O11,">$out11" or die $!;
$k = 0;
while(<I>){
	chomp;
	($a,$b) = split "\t",$_;
	unless($idd{$b}){
		$k++;
		print O11 $b."\t".$k."\n";
		$idd{$b} = $k;
	}
	if($id{$a}){
		$aa = $id{$a};
		$bb = $idd{$b};
		if($filter==1){
			if($bg{$aa}){
				$final{$aa."\t".$bb}=1;
			}
		}else{
			$final{$aa."\t".$bb}=1;
		}
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
