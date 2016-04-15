#!/usr/bin/perl -w

## default
$main_dir = $0;
@tmp = split "/",$main_dir;
$main_dir = join "/",@tmp[0..($#tmp-3)];
$src_dir  = $main_dir."/src/";
$bin_dir =  $main_dir."/bin/";

sub showhelp {
	print STDERR<<EOF;
	usage: perl src/ORA_simple_filter.pl  <enriched output file> <filtered output file> <gene set use file> [options]
		demo: perl src/ORA_simple_filter.pl demo/input_mul.txt_id.out_e-LEGO.txt demo/input_mul.txt_id.out_e-LEGO.txt_filter demo/GeneSet_human.txt_FC2_human -multi 1 
		-p_thre: p value threshold, default:0.1
		-p_col:  columns to save p value, default: 2
		-jac_thre: threshold for jaccard similarity, default: 1
		-multi: whether or not the input has multi-groups, default:0  
EOF
}	
if($#ARGV < 1){
	&showhelp();
	die;
}
###
$enrich_file   = $ARGV[0];
$output_file   = $ARGV[1];
$geneset_file  = $ARGV[2]; 
##
$multi = 0;
$p_thre = 0.1;
$jac_thre = 0.15;
$p_col = 2;
##
@cmd = @ARGV;
$cmd = join("\t",@cmd);

if($cmd =~ /-h/ || $#cmd < 1 || ($cmd !~ /-h/ && $#cmd % 2 == 1)){
	&showhelp();
	die;
}
if($#cmd > 1){
	foreach $i (3..$#cmd){
		$c = $cmd[$i];
		$recmd{$c} = $i;
	}
	foreach $i (3..$#cmd){
		$c = $cmd[$i];
		if($c eq "-multi"){$multi     = $cmd[$recmd{$c} + 1];}
		if($c eq "-p_thre"){$p_thre     = $cmd[$recmd{$c} + 1];}
		if($c eq "-p_col"){$p_col     = $cmd[$recmd{$c} + 1];}
		if($c eq "-jac_thre"){$jac_thre  = $cmd[$recmd{$c} + 1];}
	}
}
##
###########
$geneset_GSM_overlap = $geneset_file."_overlap_union.txt"; ## BIOCARTA_RELA_PATHWAY   BIOCARTA_ATM_PATHWAY    Gene_set_overlap    0.18
$geneset_module_file = $geneset_GSM_overlap."_$jac_thre\_result_iNP_MSG";
unless(-e $geneset_GSM_overlap){
	print "perl $src_dir/cal_overlap_union.pl $geneset_file\n";
	system "perl $src_dir/cal_overlap_union.pl $geneset_file";
}
unless(-e $geneset_module_file){
	$tmp = "tmp$$";
	$out_tmp = $tmp."_result_iNP_MSG";
	open ORI,"$geneset_GSM_overlap" or die $!;
	open TMP,">$tmp";
	while(<ORI>){
		chomp;
		($gs1,$gs2,$score) = split "\t";
		if($score>=$jac_thre){
			print TMP $gs1."\t".$gs2."\t".$score."\n";
		}
	}	
	close ORI;
	close TMP;
	print "$bin_dir/iNP $tmp 1 0 tmp\n";
	system "$bin_dir/iNP $tmp 1 0 tmp";
	system "mv $out_tmp $geneset_module_file";
	system "rm -rf $tmp*";
}
## read in module file
@all_gs = split "\n",`awk {'print \$2'} $geneset_file`;
$m = 0;
open M,$geneset_module_file or die $!;
while(<M>){
	chomp;
	@a = split "\t";
	if($#a < 1){
		print "Modularity:".$a[0]."\n";
		next;
	}
	$m ++;
	foreach $gs (@a){
		$belong{$gs} = $m;
	}
}
close M;
$n1 = scalar keys %belong;
print "There are $n1 gene sets in iNP results\n";
foreach $gs (@all_gs){
	if($belong{$gs}){
		next;
	}
	$m ++;
	$belong{$gs} = $m;
}
##
print "In total $m modules\n";
open E,$enrich_file or die $!;
open O,">$output_file" or die $!;
while(<E>){
	chomp;
	@a = split "\t";
	#if($. == 1){next;}
	if(/Ori_p/){next;}
	#if($a[$p_col-1] !~ /\d/){next;}
	if($multi == 1){
		$int = $a[$#a];
	}else{
		$int = "ONE";
	}
	unless($tmp{$int}){
		push(@all_int,$int);
		$tmp{$int} = 1;
	}
	$ori_gs = $a[0];
	$gs = $a[0];
	if($ori_gs =~ /(^GO:.......).*/){
		$gs = $1;
		#$des{$gs} = $ori_gs;
	}
	if($p_col == 0){
		$pv = 1;
	}else{
		$pv = $a[$p_col-1]; ######
	}
	if(abs($pv) <= $p_thre){
		$bgs = $belong{$gs};		
		if($multi == 1){
			$final{$int}{$bgs}{$ori_gs} = $pv;	
		}else{
			$final{$int}{$bgs}{$ori_gs} = $pv;	
		}
	}
}
close E;
##
open O1,">$output_file\_tab.txt" or die $!;
foreach $int (@all_int){
	foreach $bgs (keys %{$final{$int}}){
		@tmp_gs = sort{$final{$int}{$bgs}{$a} <=> $final{$int}{$bgs}{$b}}(keys %{$final{$int}{$bgs}});
		$mark_gs = $tmp_gs[0];	
		$n = "";
		foreach $gs (@tmp_gs){
			$pv = $final{$int}{$bgs}{$gs};
			$n .= $gs."_".$pv.";";
			if($multi == 1){
				print O $gs."\t".$pv."\t$bgs\t".$mark_gs."\t".$int."\n";
			}else{
				print O $gs."\t".$pv."\t$bgs\t".$mark_gs."\n";
			}
		}
		$tmp = $mark_gs."_".$final{$int}{$bgs}{$mark_gs}."\n".$n;
		$tobgs{$tmp} = $bgs;
		$cluster_out{$int}{$tmp} = $final{$int}{$bgs}{$mark_gs};
	}
}
foreach $int (@all_int){
	$k = 0;
	if($multi == 1){
		if($int ne $all_int[0]){
			print O1 "\n";
		}
		print O1 "Result for $int:\n";
	}
	foreach $tmp (sort{$cluster_out{$int}{$a} <=> $cluster_out{$int}{$b}}(keys %{$cluster_out{$int}})){
		$k++;
		print O1 "Cluster$k:".$tmp."\n";
		#print O1 "Cluster$k\_($tobgs{$tmp}):".$tmp."\n";
	}
}
close O1;
close O;
