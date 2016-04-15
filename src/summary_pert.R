#Rscript $src_dir/summary_pert.R $each $tmp_input $final_file
Args <- commandArgs(T)
file_perm <- Args[1];
out_file <- Args[2];
path <- Args[3];
##
library(multicore);
input_info <- as.numeric(system(paste(path,"/bin/extract_info ",file_perm,sep=""),intern=TRUE)[-1]);
all_gs <- as.numeric(system(paste(path,"/bin/extract_info_gs ",file_perm,sep=""),intern=TRUE)[-1]);
save(input_info,all_gs,file=out_file);
## output

