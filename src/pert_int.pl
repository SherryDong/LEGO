#!/usr/bin/perl -w
# label shuffling
use strict;
use List::Util qw/shuffle/;
if($#ARGV < 1){
	die;
}

my $a; my $b; my $w; my $i;
my %all; my %to;
my @ori; my @all; my @all_w;
my $file = $ARGV[0];
my $size = $ARGV[1];
my $per_num = $ARGV[2];
my $out_file = $ARGV[3];
# 7SK CLCF1   1
# A26C1B  1 
@ori = split "\n",`awk {'print \$1'} $file`;
map{$all{$_}=1}@ori;
@all = keys %all;
@all_w = shuffle @all;
open O,">".$out_file or die $!;
##
foreach $i (1..$per_num){
	@all_w = shuffle @all_w;
	my @tmp = @all_w[(0..($size-1))];
	my $info = join "\n",@tmp;
	foreach $a (@tmp){
		print O $a."\t$i\n";
	}
}
##
close O;
