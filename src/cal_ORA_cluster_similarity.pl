#!/usr/bin/perl -w

# this function aims to calculate similarity between cluster results

sub showhelp {
	print STDERR<<EOF;
	usage: perl src/cal_ORA_cluster_similarity.pl <input enriched clustered files,must be multi> <gene set overlap file> [options]
		demo: perl src/cal_ORA_cluster_similarity.pl demo/input_mul.txt_id.out_e-LEGO_filter_cluster.txt  
EOF
}	
if($#ARGV < 0){
	&showhelp();
	die;
}
###
$enrich_file   = $ARGV[0]; ## output of ORA_filter.pl
$geneset_over_file = $ARGV[1];
#$geneset_over_file = "demo/GeneSet_human.txt_FC2_human_overlap_union.txt"; ##
open GO,$geneset_over_file or die $!;
while(<GO>){
	chomp;
	@a = split "\t";
	$gs1 = $a[0];
	$gs2 = $a[1];
	$score{$gs1}{$gs2} = $a[2];
	$score{$gs2}{$gs1} = $a[2];
}
close GO;
##
open EF,$enrich_file or die $!;
while(<EF>){
	chomp;
	if(/^Results for (.*)/){
		$geo = $1;
		next;
	}
	if(/^Cluster(.*)\:/){
		$cluster_id = $1;
		next;
	}
	@a = split "\t";
	foreach $each (@a){
		@tmp = split "_",$each;
		$name = join "_",@tmp[(0..($#tmp-1))];
		$result{$geo}{$cluster_id}{$name} = 1;
	}
}
close EF;
##
@all_geo = keys %result;
foreach $i (@all_geo){
	@cluster_id_1 = keys %{$result{$i}};
	foreach $j (@all_geo){
		@cluster_id_2 = keys %{$result{$j}};
		$ss = 0;
		foreach $c1 (@cluster_id_1){ ## for each cluster 1, find max
			$max_score = 0;
			foreach $c2 (@cluster_id_2){
				@n1 = keys %{$result{$i}{$c1}};
				@n2 = keys %{$result{$j}{$c2}};
				## count
=cut
				undef(%tmp);map{$tmp{$_}=1}@n1;
				$D=0;
				foreach $n (@n2){
					if($tmp{$n}){
						$D++;
					}
				}
				$B = $#n1+1; $C = $#n2+1;
				#$A = $B+$C-$D;
				if($B<$C){
					$tmp_score = $D/$B;
				}else{
					$tmp_score = $D/$C;
				}
=cut				
				$t1 = 0;
				foreach $n1 (@n1){
					$tmp_max_score = 0;
					foreach $n2 (@n2){
						if($score{$n1}{$n2}){
							$tmp_max_score = ($tmp_max_score>$score{$n1}{$n2})?$tmp_max_score:$score{$n1}{$n2};
						}	
						if($n1 eq $n2){
							$tmp_max_score = 1;
						}
					}
					$t1 += $tmp_max_score;
				}
				$tmp_score = $t1/($#n1+1);
				$max_score = ($max_score>$tmp_score)?$max_score:$tmp_score;
			}
			#print $i."\t".$j."\t".$max_score."\n";
			$ss+=$max_score;
		}
		$score{$i}{$j} = $ss/($#cluster_id_1+1);
	}
}
## final score
foreach $i (@all_geo){
	foreach $j (@all_geo){
		$final = ($score{$i}{$j}+$score{$j}{$i})/2;
		print $i."\t".$j."\t".$final."\n";
	}
}


