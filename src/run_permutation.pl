#!/usr/bin/perl -w
# perl $src_dir/run_permutation.pl $interest_file\_id.out $mid_data_dir $multi
$out = $ARGV[0];
$geneset = $ARGV[1];
$net_file = $ARGV[2];
$geneset_file = $ARGV[3];
$main_dir = $ARGV[4];
$perm_times = $ARGV[5];
##
$src_dir = $main_dir."/src/";
##
$gene_id = $net_file."_id";
$mid_data_dir = $geneset."_mid_data";
##
@tmp = split "\n",`awk {'print \$9'} $out | sort | uniq`;
foreach $each (@tmp){
	if($each eq "intSize_net"){
		next;
	}
	$tmp = $each;
	if($tmp == 0){next;}
	$num{$tmp} = 1;
}
##
@num = sort{$a<=>$b}(keys %num);
print "Input size: @num\n";
##
unless(-e $mid_data_dir){
	system "mkdir $mid_data_dir";
}
foreach $each (@num){
	$final_file = "$mid_data_dir/$each.RData";
	if(-e $final_file){
		next;
	}
	$tmp_input = $geneset."_input_$each";
	print "perl $src_dir/pert_int.pl $gene_id $each $perm_times $tmp_input\n";
	system "perl $src_dir/pert_int.pl $gene_id $each $perm_times $tmp_input";
	print "perl $main_dir/LEGO_noperm.pl $net_file $geneset_file $tmp_input -multi 1 -noR 0 -min 0 -max 100000 \n";
	system "perl $main_dir/LEGO_noperm.pl $net_file $geneset_file $tmp_input -multi 1 -noR 0 -min 0 -max 100000 ";
	print "Rscript $src_dir/summary_pert.R $tmp_input\_id.out $final_file $main_dir \n";
	system "Rscript $src_dir/summary_pert.R $tmp_input\_id.out $final_file $main_dir ";
	system "rm -rf $tmp_input*";
}
