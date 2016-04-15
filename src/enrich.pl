#!usr/bin/perl -w

if($#ARGV < 2){
&showhelp();
die;
}

$database_file = $ARGV[0];
$input_file    = $ARGV[1];
$output_file   = $ARGV[2];

$item_col    = 1;
$des_col     = 2;
$des_file    = "";
$list        = 2;
$thre        = 0.01;
$background  = "no";
$min_size    = 0;
$max_size    = 10000000;
$need_adj    = 0;
$low_match   = 0;
$filter      = 0;
$filter_col  = 0;
$multi   = 0;

#### read in cmd line
@cmd = @ARGV;
$cmd = join("\t",@cmd);
#
if($cmd =~ /-h/ || $#cmd < 2 || ($cmd !~ /-h/ && $#cmd % 2 == 1)){
	&showhelp();
	die;
}
if($#cmd > 2){
	foreach $i (3..$#cmd){
		$c = $cmd[$i];
		$recmd{$c} = $i;
	}
	foreach $i (3..$#cmd){
		$c = $cmd[$i];
		if($c eq "-item"){$item_col = $cmd[$recmd{$c} + 1];}		
		if($c eq "-des"){$des_col   = $cmd[$recmd{$c} + 1];}		
		if($c eq "-des_file"){$des_file   = $cmd[$recmd{$c} + 1];}		
		if($c eq "-bg"){$background = $cmd[$recmd{$c} + 1];}		
		if($c eq "-list"){$list     = $cmd[$recmd{$c} + 1];}		
		if($c eq "-thre"){$thre     = $cmd[$recmd{$c} + 1];}		
		if($c eq "-max"){$max_size  = $cmd[$recmd{$c} + 1];}		
		if($c eq "-min"){$min_size  = $cmd[$recmd{$c} + 1];}		
		if($c eq "-adj"){$need_adj  = $cmd[$recmd{$c} + 1];}		
		if($c eq "-low"){$low_match = $cmd[$recmd{$c} + 1];}		
		if($c eq "-filter"){$filter = $cmd[$recmd{$c} + 1];}		
		if($c eq "-multi"){$multi = $cmd[$recmd{$c} + 1];}		
		if($c eq "-filter_column"){$filter_col = $cmd[$recmd{$c} + 1];}		
	}	
}

####
@method_all = ("none","BH","fdr","holm","hochberg","hommel","BH","BY");
map{$use_method_all{$_}=1}@method_all;
if($use_method_all{$need_adj}){
	$adj_method = $need_adj;
}else{
	$adj_method = $method_all[$need_adj];
}

####### read in background (if $background != 0)
if($background ne "no"){
	open B,$background or die $!;
	while(<B>){
		chomp;
		@a = split "\t";
		$exist_back{uc($a[0])}=1;
	}
	close B;
	print "Background items:";
	print scalar (keys %exist_back),"\n";
}else{
	print "no background\n";
}

###### read in des file
if($des_file ne ""){
	open D,$des_file or die $!;
	while(<D>){
		chomp;
		@a = split "\t";
		if($#a < 2){next;}
		$des{$a[0]} = $a[1]."_".$a[2];
	}
	close D;
}
###### a: whole gene numbers
open S,$database_file or die $!;
while(<S>){
	chomp($_);
	@arr  = split"\t",$_;
	$item = uc($arr[$item_col-1]);
	$des  = $arr[$des_col-1];
	if($item eq "."){next;}
	if($des_file ne ""){
		if($des{$des}){$des = $des."_".$des{$des};}
	}
	if($filter){
		$use_filter = uc($arr[$filter_col-1]);
		if($filter ne $use_filter){
			next;
		}
	}
	if($des !~ /\S/){next;}
	if($background ne "no"){
		if(exists($exist_back{$item})){
			$exist{$item}=1;
			if($des !~ /;/){	
				$all{$des}{$item} = 1;	
			}else{
				@des = split ";",$des;
				map{$all{$_}{$item} = 1}@des;	
			}
		}
		if($des !~ /;/){	
			$all_size{$des}{$item} = 1;	
		}else{
			@des = split ";",$des;
			map{$all_size{$_}{$item} = 1}@des;	
		}
	}
	if($background eq "no"){
		$exist{$item}=1;
		if($des !~ /;/){	
			$all{$des}{$item} = 1;	
		}else{
			@des = split ";",$des;
			map{$all{$_}{$item} = 1}@des;	
		}
	}
}
close S;
$a = scalar (keys %exist);
print "All items for analysis:".$a,"\n";
if($background eq "no"){
	%all_size = %all;
}
foreach $each_des (keys %all_size){
	$b = scalar (keys %{$all_size{$each_des}});
	if($b < $min_size || $b > $max_size){
		$no{$each_des} = 1;
	}
}
#######input list c--> whole gene numbers in the list
open IN,$input_file or die $!;
open OUT,">$output_file" or die $!;
while(<IN>){
	chomp($_);
	@arr = split"\t",$_;
	$gene = uc($arr[0]);
	if($multi){
		$group = $arr[1];
		if($background ne "no"){
		if($exist_back{$gene}){
			$all_count{$group}{$gene} = 1;
		}}else{
			$all_count{$group}{$gene} = 1;
		}
		$mark{$group}{$gene} = 1;
	}else{
		$mark{1}{$gene}=1;
		if($background ne "no"){
			if($exist_back{$gene}){
				$all_count{$gene} = 1;
			}
		}else{
			$all_count{$gene} = 1;
		}
	}
}
close IN;
##
$ng = scalar (keys %mark);
print "Number of groups:".$ng,"\n";

#########enrichment --> b-> gene numbers in a GO / d-> gene numbers in a GO in the list
if($multi){
	if($list == 2){
		print OUT "##Name\tTotal_Item\tNum_item\tNum_list\tNum_list_item\tOri_p\tAdj_p\tEnrich_score\tGroup_ID\n";
	}else{
		print OUT "##Name\tTotal_Item\tNum_item\tNum_list\tNum_list_item\tOri_p\tAdj_p\tEnrich_score\tGene_list\tGroup_ID\n";
	}
}else{
	if($list == 2){
		print OUT "##Name\tTotal_Item\tNum_item\tNum_list\tNum_list_item\tOri_p\tAdj_p\tEnrich_score\n";
	}else{
		print OUT "##Name\tTotal_Item\tNum_item\tNum_list\tNum_list_item\tOri_p\tAdj_p\tEnrich_score\tGene_list\n";
	}
}
##
$A = $a;
if($background ne "no"){
	$A = scalar keys %exist_back;
	unless($multi){$C = scalar keys %all_count;}
}
foreach $group (keys %mark){ ## for each group
	undef(%tofour);
	undef(%toinfo);
	undef(%togene);
	undef(%result);
	$tmp = "tmp".$$;
	open TR,">$tmp";
	$i = 0;
	$j = 0;
	foreach $each_des (keys %all){
		if($no{$each_des}){next;}
		$b = scalar (keys %{$all{$each_des}});
		undef(@gene_all);	
		@item_all = keys %{$all{$each_des}};
		$c = scalar keys %{$mark{$group}};
		if($background ne "no"){
			if($multi){
				$c = keys %{$all_count{$group}};
			}else{
				$c = $C;
			}
		}
		@gene_all = grep{$mark{$group}{$_}}@item_all;
		$d = $#gene_all + 1;
		if($d < $low_match){
			next;
		}
		##
		$a = $A;
		@input = ($a,$b,$c,$d);
		$four = join("\t",@input);
		#print $four."\n";
		$tofour{$i} = $four;
		$toinfo{$i} = $each_des."\t".$four;
		$togene{$i} = join(",",@gene_all);
		$i ++;
		print TR $four."\n";
		#if($a*$d <= $b*$c){
		#	next;
		#}
		$j ++;
	}

	if($j == 0){
		print "No enriched ! Change the parameters & try!\n";
		unlink ($tmp) if (-f $tmp); 
		next;
	}
################## output result

	@all_out = fisher_test();

	foreach $i (0..$#all_out){
		$OUT  = $all_out[$i];
		$info = $toinfo{$i};
		$gene_all = $togene{$i};
		$four = $tofour{$i};
		($a,$b,$c,$d) = split"\t",$four;
		#if($a*$d <= $b*$c){
		#	next;
		#}
		$ratio = ($d/$c)/($b/$a);
		$out = (split"_",$OUT)[0];
		$adj = (split"_",$OUT)[1];
		if($out eq "null"){next;}
		if($adj > $thre){next;}
		$p   = $out;
		$out = $out*(1+$i/10000000+1/($$));
		if($list == 2){
			$result{$out}=join("\t",$info,$p,$adj,$ratio);
		}else{
			$result{$out}=join("\t",$info,$p,$adj,$ratio,$gene_all);
		}
	}

	foreach $each(sort {$a <=> $b}(keys %result)){
		if($multi){
			print OUT $result{$each}."\t$group\n";
		}else{
			print OUT $result{$each}."\n";
		}
	}
}
close OUT;

#########
sub fisher_test {
	my ($each, $p);
	my $out = "$tmp.R";
	open TMR,">$out" or die $!;

	print TMR<<CMD;
	all <- read.table(file = "$tmp",sep = "\t")
		all <- as.matrix(all)	
		p   <- apply(all,1,function(x){
				fisher.test(matrix(x, nrow=2),alternative="greater")\$p.value
				})
	adj.p <- p.adjust(p,method="$adj_method")
		result <- cbind(p,adj.p)
		print(result,quote=F,col.names=F,row.names=F)
CMD
		close (TMR);
	my $txt      = `R --no-save < $out`;
	my @tmp1     = split">",$txt;
	my $use_info = $tmp1[6];
	my @tmp2     = split"\n",$use_info;
	my $i;
	my @out;
	my @arr;
	foreach $i (2..$#tmp2){
		@arr = split" ",$tmp2[$i];
		$each = $arr[1]."_".$arr[2];
		push(@out,$each);
	}
	unlink ($out) if (-f $out); 
	unlink ($tmp) if (-f $tmp); 
	return @out; 
}

####
sub showhelp {
	print STDERR<<EOF;
	enrich.pl  (v1.1)

		Usage:  enrich.pl  <database file> <input_file_name> <output_file_name>  [options]

		database file :  database file for enrichment analyais 

		input_file_name:   input file name

		output_file_name: 	output file name

		Options:  

		-h  Show help 

		-item   (item column number,default = 1)

		-des    (des column number, default = 2) 

		-des_file (files to save description for genesets)

		-bg     (bacground_file_name,default=no background file)

		-list   (1-y|2-n,default=n)

		-thre   (threshold,default = 0.01)

		-max    (maximun size for des, default = no limitation)

		-min    (minmun size for des, default = no limitation)

		-adj    (need adj p value or not ,0-no,1-BH,2-FDR,default = no)

		-low    (low existence for match,default = 0)
		 
		-multi  (whether the input has multi groups, default=0)

		-filter ( use what filter or not, default = 0 (do not use))

		-filter_column (filter column)
EOF
}

