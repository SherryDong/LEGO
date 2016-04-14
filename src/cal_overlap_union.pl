#!/usr/bin/perl -w

$geneset_use_file = $ARGV[0];
$output_file = $geneset_use_file."_overlap_union.txt";
open O,">".$output_file or die $!;
# BIOCARTA_RELA_PATHWAY   BIOCARTA_AKT_PATHWAY    Gene_set_overlap    0.37
open GS,$geneset_use_file or die $!;
while(<GS>){
	chomp;
	($gene,$gs) = split "\t";			
	$all{$gs}{$gene} = 1;
	#$total{$gene} = 1;
}
close GS;
#$A = keys %total;
@all_gs = keys %all;
$gs_num = $#all_gs+1;
$total_adj = $gs_num*($gs_num-1)/2;
print $total_adj."\n";
foreach $a (0..($#all_gs-1)){
	$gs1 = $all_gs[$a];
	$B = keys %{$all{$gs1}};
	@all_gene = keys %{$all{$gs1}};
	foreach $b (($a+1)..$#all_gs){
		$gs2 = $all_gs[$b];
		$C = keys %{$all{$gs2}};
		$D = 0;
		foreach $g (@all_gene){
			if($all{$gs2}{$g}){
				$D ++;
			}
		}
		if($D == 0){next;}
		 $A = $B+$C-$D;
       		 $jac = $D/$A;
       		 print O $gs1."\t".$gs2."\t$jac\n";
	}
}
close O;
###

