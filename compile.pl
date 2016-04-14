#! /usr/bin/perl -w 


$main_dir = `pwd`;
$src_dir = $main_dir."/src/"; 
$src_dir =~s/\s//g;
$exec_dir = $src_dir; 
$exec_dir =~s/src/bin/; 
if (!-d $exec_dir) { 
	mkdir (0755, $exec_dir); 
}
## LEGO
system "gcc -o $exec_dir/LEGO_pre $src_dir/LEGO_pre.c -lm";
system "gcc -o $exec_dir/LEGO $src_dir/LEGO.c -lm";
system "gcc -o $exec_dir/LEGO_mul $src_dir/LEGO_mul.c -lm";
system "gcc -o $exec_dir/LEGO_bg $src_dir/LEGO_bg.c -lm";
system "gcc -o $exec_dir/LEGO_mul_bg $src_dir/LEGO_mul_bg.c -lm";
system "g++ -o $exec_dir/extract_info $src_dir/extract_info.cpp -lm";
system "g++ -o $exec_dir/extract_info_gs $src_dir/extract_info_gs.cpp -lm";
##

$msg_dir = $src_dir."/MSG/"; 

$inp = "$src_dir/iNP.c"; 
$new_code = $inp.".new"; 
open (N, ">$new_code") || die; 
open (C, $inp) || die; 
while (<C>) { 
	if (/(#define)\s(BIN_DIR)\s(\S+)/) {
		print N $1, " ", $2, " \"", $exec_dir, "/\"\n"; 
	} else { 
		print N $_; 
	}
}
close (C); 
close (N); 
rename($new_code, $inp); 
system "gcc -o $exec_dir/iNP $src_dir/iNP.c -lm"; 

chdir ($msg_dir); 
$code1 = "$src_dir/MSG.cpp"; 
$code2 = "$src_dir/VM.cpp"; 
$new_code = $code1.".new"; 
open (N, ">$new_code") || die; 
open (C, $code1) || die; 
while (<C>) { 
	if (/(#define)\s(VM_dir)\s(\S+)/) {
		print N $1, " ", $2, " \"", $exec_dir, "/\"\n"; 
	} else { 
		print N $_; 
	}
}
close (C); 
close (N); 
rename($new_code, $code1); 
system "g++ $code1 -o MSG.out"; 
system "g++ $code2 -o VM.out"; 
system "mv MSG.out $exec_dir"; 
system "mv VM.out $exec_dir"; 

