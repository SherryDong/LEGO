#!/usr/bin/perl -w

use strict;
use Cwd;

sub showhelp {
	print STDERR<<EOF;
	perl LEGO_noperm.pl <network file> <geneset file> <interest file> [options]
	
	Options:  

	-h/-help: show help message.
	-pre: <whether or not to pre-run to generate mid files,1-yes,0-no(e.g:0)> 
	-multi: <multi or not:1-yes,0-no(e.g:0)> 
	-bgNE: <bgNE for the network,e.g:0.25> 
	-min: <minimun gene set size,e.g:5> 
	-max: <maximum gene set size,e.g:10000> 
	-adj: <adjusted methods,e.g:fdr> 
	-p_thre: <p value cutoff,e.g:0.05> 
	-fisher: <whether or not run fisher exact test, 1-yes,0-no,e.g: 0> 
	-filter: <whether or not run result cluster and filter step, 1-yes, 0-no, e.g: 0> 
	-bg_file: <background file, e.g: bg.txt, if do not have background file, leave it blank> \n";
	-noR: <whether or not to run the following R script> 1: yes, 0, no; default: 1;

#############################################################################################
EOF
}
#
if($#ARGV < 2){
	&showhelp();
	die;
}
############################################## default paras

my $network_file  = $ARGV[0]; ## network file
my $geneset_file  = $ARGV[1]; ## gene set file
my $interest_file = $ARGV[2]; ## interesting file
##
my $cmd;
my @cmd;
my $network_mark;
my $geneset_use_file;
my $c;
my $i;
my %recmd;
my @tmp;
my $main_dir = $0;
@tmp = split "/",$main_dir;
if($#tmp<1){
	$main_dir = getcwd()."/".$0;
	@tmp = split "/",$main_dir;
}
$main_dir = join "/",@tmp[0..($#tmp-1)];
my $src_dir = "$main_dir/src/";
my $exec_dir = "$main_dir/bin/";
## options
my $pre_run = 0;
my $min=5;
my $max=10000;
my $multi = 0;
my $bgNE = 0.25;
my $adj_meth = "fdr";
my $p_thre = 0.05;
my $fisher = 0;
my $filter = 0;
my $bg_file = "";
my $no_R = 1;

###################### 
@cmd = @ARGV;
$cmd = join("\t",@cmd);

if($cmd =~ /-h/ || $#cmd < 1 || ($cmd !~ /-h/ && $#cmd % 2 == 1)){
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
		if($c eq "-pre"){$pre_run  = $cmd[$recmd{$c} + 1];}
		if($c eq "-min"){$min  = $cmd[$recmd{$c} + 1];}
		if($c eq "-max"){$max  = $cmd[$recmd{$c} + 1];}
		if($c eq "-multi"){$multi  = $cmd[$recmd{$c} + 1];}
		if($c eq "-bgNE"){$bgNE  = $cmd[$recmd{$c} + 1];}
		if($c eq "-adj"){$adj_meth  = $cmd[$recmd{$c} + 1];}
		if($c eq "-p_thre"){$p_thre  = $cmd[$recmd{$c} + 1];}
		if($c eq "-fisher"){$fisher  = $cmd[$recmd{$c} + 1];}
		if($c eq "-filter"){$filter  = $cmd[$recmd{$c} + 1];}
		if($c eq "-bg_file"){$bg_file  = $cmd[$recmd{$c} + 1];}
		if($c eq "-noR"){$no_R  = $cmd[$recmd{$c} + 1];}
 }
}

# run program
my $date = `date`;
chomp($date);
print "\n###### Run begin:$date ######\n";
@tmp = split "\/",$network_file;
$network_mark = $tmp[$#tmp];
$network_mark =~ s/\.txt//g;
$geneset_use_file = $geneset_file."_".$network_mark;
print "\n*Run pre-programs to generate mid-files:\n";
if($pre_run){
	system "perl $src_dir/parse_net.pl $network_file";
	system "perl $src_dir/parse_gs.pl $geneset_file $geneset_use_file $network_file\_id";
	print "\n*Run LEGO_pre to generate mid files:\n\t$exec_dir/LEGO_pre $network_file\_id_net $geneset_use_file\_id_gs $bgNE\n";
	system "$exec_dir/LEGO_pre $network_file\_id_net $geneset_use_file\_id_gs $bgNE";
}else{
	unless(-e "$network_file\_id"){
		system "perl $src_dir/parse_net.pl $network_file";
		system "perl $src_dir/parse_gs.pl $geneset_file $geneset_use_file $network_file\_id";
		print "\n*Run LEGO_pre to generate mid files:\n\t$exec_dir/LEGO_pre $network_file\_id_net $geneset_use_file\_id_gs $bgNE\n";
		system "$exec_dir/LEGO_pre $network_file\_id_net $geneset_use_file\_id_gs $bgNE";
	}else{ 
		unless(-e "$geneset_use_file\_id_gs"){system "perl $src_dir/parse_gs.pl $geneset_use_file $network_file\_id";}
		unless(-e "$geneset_use_file\_id_gs_NW"){system "$exec_dir/LEGO_pre $network_file\_id_net $geneset_use_file\_id_gs $bgNE";}
	}
}
## no bg
if($bg_file){
	if($multi){
		system "perl $src_dir/parse_int_multi_bg.pl $interest_file $bg_file $geneset_use_file\_gene_id 1";
		$cmd = "$exec_dir/LEGO_mul_bg $geneset_use_file\_id_gs $geneset_use_file\_id_gs_NW $geneset_use_file\_id_gs_CS $geneset_use_file\_id_gs_GS $interest_file\_id $min $max $bg_file\_id";
		print "\n*Run main program:\n\t".$cmd."\n";
		system $cmd;
	}else{
		system "perl $src_dir/parse_int_bg.pl $interest_file $bg_file $geneset_use_file\_gene_id 1";
		$cmd = "$exec_dir/LEGO_bg $geneset_use_file\_id_gs $geneset_use_file\_id_gs_NW $geneset_use_file\_id_gs_CS $geneset_use_file\_id_gs_GS $interest_file\_id $min $max $bg_file\_id";
		print "\n*Run main program:\n\t".$cmd."\n";
		system $cmd;
	}
	if($fisher){
		print "\n*Run Fisher's Exact test:\n\tperl $src_dir/enrich.pl $geneset_use_file $interest_file $interest_file\_fisher.txt -max $max -min $min -list 1 -thre $p_thre -adj $adj_meth -multi $multi -bg $bg_file\n";
		system "perl $src_dir/enrich.pl $geneset_use_file $interest_file $interest_file\_fisher.txt -max $max -min $min -list 1 -thre $p_thre -adj $adj_meth -multi $multi -bg $bg_file";
	}
}else{
	if($multi){
		system "perl $src_dir/parse_int_multi.pl $interest_file $geneset_use_file\_gene_id";
		$cmd = "$exec_dir/LEGO_mul $geneset_use_file\_id_gs $geneset_use_file\_id_gs_NW $geneset_use_file\_id_gs_CS $geneset_use_file\_id_gs_GS $interest_file\_id $min $max";
		print "\n*Run main program:\n\t".$cmd."\n";
		system $cmd;
	}else{
		system "perl $src_dir/parse_int.pl $interest_file $geneset_use_file\_gene_id";
		$cmd = "$exec_dir/LEGO $geneset_use_file\_id_gs $geneset_use_file\_id_gs_NW $geneset_use_file\_id_gs_CS $geneset_use_file\_id_gs_GS $interest_file\_id $min $max";
		print "\n*Run main program:\n\t".$cmd."\n";
		system $cmd;
	}
	if($fisher){
		print "\n*Run Fisher's Exact test:\n\tperl $src_dir/enrich.pl $geneset_use_file $interest_file $interest_file\_fisher.txt -max $max -min $min -list 1 -thre $p_thre -adj $adj_meth -multi $multi\n";
		system "perl $src_dir/enrich.pl $geneset_use_file $interest_file $interest_file\_fisher.txt -max $max -min $min -list 1 -thre $p_thre -adj $adj_meth -multi $multi";
	}
}
if($no_R==0){
	die;
}
system "perl $src_dir/interpret.pl $interest_file\_id.out  $geneset_use_file\_id $adj_meth $p_thre $multi";
### filter
if($filter){
	unless(-e "$geneset_use_file\_overlap.txt"){
		print "\n*Run cal_overlap to calculate overlap value between gene set pairs:\n\tperl $src_dir/cal_overlap_noperm.pl $geneset_use_file";
		system "perl $src_dir/cal_overlap.pl $geneset_use_file";
	}
	if($bg_file){
		print "\n*Run ORA_filter to get filtered and clustered enriched gene sets:\n\tperl $src_dir/ORA_filter.pl $interest_file $interest_file\_id.out_e-LEGO.txt $interest_file\_id.out_e-LEGO_filter $geneset_use_file $network_file\_id_net -multi $multi -bg $bg_file \n";	
		system "perl $src_dir/ORA_filter.pl $interest_file $interest_file\_id.out_e-LEGO.txt $interest_file\_id.out_e-LEGO_filter $geneset_use_file $network_file\_id_net -multi $multi -bg $bg_file ";	
	}else{
		print "\n*Run ORA_filter to get filtered and clustered enriched gene sets:\n\tperl $src_dir/ORA_filter.pl $interest_file $interest_file\_id.out_e-LEGO.txt $interest_file\_id.out_e-LEGO_filter $geneset_use_file $network_file\_id_net -multi $multi \n";	
		system "perl $src_dir/ORA_filter.pl $interest_file $interest_file\_id.out_e-LEGO.txt $interest_file\_id.out_e-LEGO_filter $geneset_use_file $network_file\_id_net -multi $multi";	
	}
}
### finish
$date = `date`;
chomp($date);
print "###### Run end:$date ######\n\n";

