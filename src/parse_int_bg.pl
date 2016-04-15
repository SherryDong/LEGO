#!/usr/bin/perl -w
# This script read in interesting gene list file and convert to integer ID

if($#ARGV < 1){
	print "usage: perl src/parse_int_bg.pl <interesting gene list file> <background file> <gene Id file> <use bg filter or not> \n";
	die;
}
$input = $ARGV[0]; ## input gene set file
$bg    = $ARGV[1]; ## input gene set file
$id_f  = $ARGV[2]; ## input gene Id file
$filter = $ARGV[3];
$out1  = $input."_id";
$out2  = $bg."_id";
##
open ID,$id_f or die $!;
while(<ID>){
	chomp;
	($g,$i) = split "\t";
	$id{$g} = $i;
}
close ID;
## read in bg file
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
## read in input file
open I,$input or die $!;
open O1,">$out1" or die $!;
while(<I>){
	chomp;
	@a = split "\t";
	$a = $a[0];
	if($id{$a}){
		$aa = $id{$a};
		if($filter==1){
			if($bg{$aa}){
				$final{$aa} = 1;
			}
		}else{
			$final{$aa} = 1;
		}
	}else{
		$no{$a}=1;
	}
}
close I;
@n = sort{$a<=>$b}(keys %final);
$n = join "\n",@n;
print O1 $n."\n";
close O1;
@no = keys %no;$no=$#no+1;
if($no>0){
	print "There are $no genes from interesting gene list or background file not found in the $id_f!\n";
}
