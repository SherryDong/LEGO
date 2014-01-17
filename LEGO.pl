#!/usr/bin/perl -w

use strict;

if($#ARGV < 2){
	print "usage: perl LEGO.pl <network file> <geneset file> <interest file> <multi or not:1-yes,0-no(e.g:0)> <bgNE,e.g:0.25> <minimun gene set size,e.g:5> <maximum gene set size,e.g:10000> <adjusted methods,e.g:fdr> <p value cutoff,e.g:0.05>\n";
	die;
}
#
my $compile = 0;
my $pre_run = 0;
# input
my $min=5;
my $max=10000;
my $cmd;
my $network_file  = $ARGV[0];
my $geneset_use_file  = $ARGV[1];
my $interest_file = $ARGV[2];
my $multi = 0;
if($#ARGV > 2){$multi = $ARGV[3];}
my $bgNE = 0.25;
if($#ARGV > 3){$bgNE = $ARGV[4];}
if($#ARGV > 4){$min = $ARGV[5];}
if($#ARGV > 5){$max = $ARGV[6];}
my $adj_meth = "fdr";
my $p_thre = 0.05;
if($#ARGV > 6){$adj_meth = $ARGV[7];}
if($#ARGV > 7){$p_thre = $ARGV[8];}
# compile 
if($compile){
	my $main_dir = `pwd`; 
	$main_dir =~s/\s//g;
	my $src_dir  = $main_dir."/src";
	my $bin_dir = $main_dir."/bin"; 
	if (!-d $bin_dir) { 
		system "mkdir $bin_dir"; 
	}
	system "gcc -o $bin_dir/LEGO $src_dir/LEGO.c -lm";
	system "gcc -o $bin_dir/LEGO_mul $src_dir/LEGO_mul.c -lm";
	system "gcc -o $bin_dir/LEGO_pre $src_dir/LEGO_pre.c -lm";
}
# run program
my $date = `date`;
print "\nRun begin:$date";
if($pre_run){
	system "perl src/parse_net.pl $network_file";
	system "perl src/parse_gs.pl $geneset_use_file $network_file\_id";
	system "./bin/LEGO_pre $network_file\_id_net $geneset_use_file\_id_gs $bgNE";
}else{
	unless(-e "$network_file\_id"){system "perl src/parse_net.pl $network_file";} 
	unless(-e "$geneset_use_file\_id_gs"){system "perl src/parse_gs.pl $geneset_use_file $network_file\_id";}
	unless(-e "$geneset_use_file\_id_gs_NW"){system "./bin/LEGO_pre $network_file\_id_net $geneset_use_file\_id_gs $bgNE";}
}
if($multi){
	system "perl src/parse_int_multi.pl $interest_file $geneset_use_file\_gene_id";
	$cmd = "./bin/LEGO_mul $geneset_use_file\_id_gs $geneset_use_file\_id_gs_NW $geneset_use_file\_id_gs_CS $geneset_use_file\_id_gs_GS $interest_file\_id $min $max";
	print $cmd."\n";
	system $cmd;
}else{
	system "perl src/parse_int.pl $interest_file $geneset_use_file\_gene_id";
	$cmd = "./bin/LEGO $geneset_use_file\_id_gs $geneset_use_file\_id_gs_NW $geneset_use_file\_id_gs_CS $geneset_use_file\_id_gs_GS $interest_file\_id $min $max";
	print $cmd."\n";
	system $cmd;
}
system "perl src/interpret.pl $interest_file\_id.out  $geneset_use_file\_id $adj_meth $p_thre $multi";
$date = `date`;
print "Run end:$date\n";
