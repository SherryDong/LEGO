#!/usr/bin/perl -w
open D,"/home/new/HY/Leo/LEGO_server/data/GeneSet/GeneSet_des.txt" or die $!;
while(<D>){
	chomp;
	@a = split "\t";
	if($#a < 2){next;}
	$des{$a[0]} = $a[1]."_".$a[2];
}
close D;
open I,$ARGV[0];
while(<I>){
	chomp;
	@a = split "\t";
	#$pv = $a[3];
	#if($pv>0.05){next;}
	$gs = $a[0];
	# Cluster419_(265):GO:0042634_0.095617_
	if($gs =~ /(GO:.......)/){
		$gs = $gs."_".$des{$gs};
		$a[0] = $gs;
	}
	$n = join "\t",@a;
	print $n."\n";
}
close I;
