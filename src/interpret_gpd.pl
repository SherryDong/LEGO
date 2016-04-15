#!/usr/bin/perl -w
# This script read in output file, convert to original name
# usage: perl src/interpret.pl <out file> <gene set ID file>
# output: 
# <out file>_txt : e.g : test_data/Int_Gene_list.txt_id.out.txt

use Cwd;
$main_dir = $0;
@tmp = split "/",$main_dir;
$main_dir = join "/",@tmp[0..($#tmp-1)];

if($#ARGV < 1){
	print "usage: perl src/interpret.pl <out file> <gene set ID file> <adjusted method> <p-value threshold> <multi> \n";
	die;
}
$input_id = $ARGV[0]; ## input out file
$each_gs  = $ARGV[1]; ## input gene set Id file
$id_s  = $each_gs."_id";
$gene_id = $each_gs."_gene_id";
$input = $input_id.".out";
$out1  = $input."_LEGO.txt"; ##
$met   = $ARGV[2];
$thre  = $ARGV[3];
$multi = $ARGV[4];
$mid_prefix = $ARGV[5];
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
$gs_id_num = 0;
open ID,$id_s or die $!;
while(<ID>){
	chomp;
	($g,$i) = split "\t";
	$ids{$i} = $g;
	$gs_id_num ++;
}
close ID;
##
## convert Z to P and do adjustment
open D,"/home/new/HY/Leo/LEGO_server/data/GeneSet/GeneSet_des.txt" or die $!;
while(<D>){
	chomp;
	@a = split "\t";
	if($#a < 2){next;}
	$des{$a[0]} = $a[1]."_".$a[2];
}
close D;
##
$tmp = "$input\_tmp".$$;
print "Rscript $main_dir/pval_gpd.R $input $thre $met $tmp $multi $mid_prefix\n";
system "Rscript $main_dir/pval_gpd.R $input $thre $met $tmp $multi $mid_prefix";
########################## transfer id to gene
@use_input = split "\n",`cat $input_id`; ## input gene id
map{$use_input{$_}=1}@use_input; ## use input gene id 
$use_input = join "\t",@use_input;
$use_input_num = $#use_input + 1;
$gene_id_num = 0;
open GI,$gene_id or die $!;
while(<GI>){
	chomp;
	$gene_id_num++;
	@a = split "\t";
	$id2gene{$a[1]} = $a[0];
}
close GI;
open GS,$each_gs."_id_gs" or die $!;
while(<GS>){
	chomp;
	@a = split "\t";
	$gs = $a[1]; $gene = $a[0];
	$belong{$gs}{$gene} = 1;
}
close GS;
######################### get results
$tmp_out = $tmp."_EdgeResults.txt";
@ori_gs = split "\n",`awk {'print \$1'} $tmp_out`;
$use_gs_num = $#ori_gs+1;
$use_gs = join "\t",@ori_gs;
##
$txt = `$main_dir/extract_CS $each_gs\_id_gs_CS $gs_id_num $use_gs_num $use_gs $gene_id_num $use_input_num $use_input`; ## neighbor: gsid\tgeneid\tweight
@txt = split "\n",$txt;
undef(%score_ov);
undef(%score_nb);
foreach $each (@txt){
	($ori_gs,$geneid,$weight) = split "\t",$each;
	if($belong{$ori_gs}{$geneid}){
		$score_ov{$ori_gs}{$id2gene{$geneid}} = $weight;
	}else{
		$score_nb{$ori_gs}{$id2gene{$geneid}} = $weight;
	}
}
foreach $each_gs (@ori_gs){
	@ov_result = sort{$score_ov{$each_gs}{$b}<=>$score_ov{$each_gs}{$a}}(keys %{$score_ov{$each_gs}});
	@nb_result = sort{$score_nb{$each_gs}{$b}<=>$score_nb{$each_gs}{$a}}(keys %{$score_nb{$each_gs}});
	$ov_result{$each_gs} = join ",",@ov_result;
	$nb_result{$each_gs} = join ",",@nb_result;
}

#########################
open I1,$tmp_out or die $!;
open O1,">$out1" or die $!;
while(<I1>){
	chomp;
	if($. == 1){
		print O1 $_."\tEssential overlapped genes (sorted by weight)\tEssential neighbor genes (sorted by weight)\n";
		next;
	}
	@a = split "\t";
	@v = @a[(1..$#a)];
	if($multi > 0){
		$v[$#v] = $toID{$v[$#v]};	
	}
	$v = join "\t",@v;
	$ori_gs = $a[0];
	$gs = $ids{$ori_gs};
	if($des{$gs}){$gs=$gs."_".$des{$gs};}
	print O1 $gs."\t".$v."\t$ov_result{$ori_gs}\t$nb_result{$ori_gs}\n";
}
close I1;
close O1;
##
print "Finish! Check output file: $out1\n";
#print "Finish! Check output file: $out1\n\t\t\t$out2\n";
system "rm -rf $tmp";
system "rm -rf $tmp\_EdgeResults.txt";
