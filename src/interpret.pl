#!/usr/bin/perl -w
# This script read in output file, convert to original name
# usage: perl src/interpret.pl <out file> <gene set ID file>
# output: 
# <out file>_txt : e.g : test_data/Int_Gene_list.txt_id.out.txt
use Cwd;

if($#ARGV < 1){
	print "usage: perl src/interpret.pl <out file> <gene set ID file> <adjusted method> <p-value threshold> <multi> \n";
	die;
}
$main_dir = $0;
@tmp = split "/",$main_dir;
$main_dir = join "/",@tmp[0..($#tmp-1)];

$input = $ARGV[0]; ## input out file
$id_s  = $ARGV[1]; ## input gene set Id file
$out1  = $input."_LEGO.txt";
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
## convert Z to P and do adjustment
$tmp = "$input\_tmp".$$;
print "Rscript $main_dir/pval.R $input $thre $met $tmp $multi \n";
system "Rscript $main_dir/pval.R $input $thre $met $tmp $multi ";
open I1,$tmp."_EdgeResults.txt" or die $!;
open O1,">$out1" or die $!;
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
	print O1 $gs."\t".$v."\n";
}
close I1;
print "Finish! Check output file: $out1\n";
#print "Finish! Check output file: $out1\n\t\t\t$out2\n";
system "rm -rf $tmp";
system "rm -rf $tmp\_EdgeResults.txt";
