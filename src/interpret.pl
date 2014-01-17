#!/usr/bin/perl -w
# This script read in output file, convert to original name
# usage: perl src/interpret.pl <out file> <gene set ID file>
# output: 
# <out file>_final : e.g : test_data/Int_Gene_list.txt_id.out.final

if($#ARGV < 1){
	print "usage: perl src/interpret.pl <out file> <gene set ID file> <adjusted method> <p-value threshold> <multi> \n";
	die;
}
$main_dir = "/home/new/HY/Leo/LEGO_server/";
$input = $ARGV[0]; ## input out file
$id_s  = $ARGV[1]; ## input gene set Id file
$out1  = $input."_n-LEGO.final";
$out2  = $input."_e-LEGO.final";
$met   = $ARGV[2];
$thre  = $ARGV[3];
$multi = $ARGV[4];
if($multi>0){
	$int_id = $input;
	$int_id =~ s/_id.*/_id_2/g;
	open II,$int_id or die $!;
	while(<II>){
		chomp;
		($a,$b) = split "\t";
		$toID{$b} = $a;
	}
	close II;
}
##
open ID,$id_s or die $!;
while(<ID>){
	chomp;
	($g,$i) = split "\t";
	$ids{$i} = $g;
}
close ID;
##
open D,$main_dir."/data/GeneSet/GeneSet_des.txt" or die $!;
while(<D>){
	chomp;
	@a = split "\t";
	if($#a < 2){next;}
	$des{$a[0]} = $a[1]."_".$a[2];
}
close D;
## convert Z to P and do adjustment
$tmp = "$input\_tmp".$$;
print "Rscript $main_dir/src/pval.R $input $thre $met $tmp $multi\n";
system "Rscript $main_dir/src/pval.R $input $thre $met $tmp $multi";
open I1,$tmp."_NodeResults.txt" or die $!;
# "GOID","Z_score","Ori_p","Adj_p","OverlapSize"
open I2,$tmp."_EdgeResults.txt" or die $!;
open O1,">$out1" or die $!;
open O2,">$out2" or die $!;
while(<I1>){
	chomp;
	if($. == 1){
		print O1 $_."\n";
		next;
	}
	@a = split "\t";
	@v = @a[(1..$#a)];
	if($multi > 0){
		$v[$#v] = $toID{$v[$#v]};	
	}
	$v = join "\t",@v;
	$gs = $ids{$a[0]};
	if($des{$gs}){
		$des = $des{$gs};
		print O1 $gs."_".$des."\t".$v."\n";
	}elsif($gs =~ /(GO:.......)/){
		$des = $des{$1};
		if($des){
			print O1 $gs."_".$des."\t".$v."\n";
		}else{
			print O1 $gs."\t".$v."\n";
		}
	}else{
		print O1 $gs."\t".$v."\n";
	}
}
close I1;
close O1;
while(<I2>){
	chomp;
	if($. == 1){
		print O2 $_."\n";
		next;
	}
	@a = split "\t";
	@v = @a[(1..$#a)];
	if($multi > 0){
		$v[$#v] = $toID{$v[$#v]};	
	}
	$v = join "\t",@v;
	$gs = $ids{$a[0]};
	if($des{$gs}){
		$des = $des{$gs};
		print O2 $gs."_".$des."\t".$v."\n";
	}elsif($gs =~ /(GO:.......)/){
		$des = $des{$1};
		if($des){
			print O2 $gs."_".$des."\t".$v."\n";
		}else{
			print O2 $gs."\t".$v."\n";
		}
	}else{
		print O2 $gs."\t".$v."\n";
	}
}
close I2;close O2;
print "Finish! Check output file: $out1\n\t\t\t$out2\n";
system "rm -rf $tmp";
system "rm -rf $tmp\_NodeResults.txt";
system "rm -rf $tmp\_EdgeResults.txt";
