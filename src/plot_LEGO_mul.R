#!/home/sherry/bin/Rscript

Args <- commandArgs(T)
if(length(Args)<6){
	print("Usage: plot_LEGO.R <input interesting gene list file> <enriched output file> <output file prefix> <network file name> <network+gene set prefix,e.g: demo/GeneSet_human.txt_FC2_human> <plot gene set name, can be comma saperated> \n");
	quit("no");
}
##
int_file    <- Args[1]
input_file  <- Args[2] ## 
output_file <- Args[3]
net_file    <- Args[4]
mark_prefix <- Args[5]
gs_for_plot <- Args[6]
gs_for_plot <- unlist(strsplit(gs_for_plot,",")); 
###for example:
#int_file <- "demo/input.txt";
#input_file <- "demo/input.txt_id.out_e-LEGO.txt";
#output_file <- "demo/test_plot";
#net_file <- "demo/FC2_human.txt";
#mark_prefix <- "demo/GeneSet_human.txt_FC2_human";
#gs_for_plot <- c("REACTOME_G1_S_SPECIFIC_TRANSCRIPTION","REACTOME_KINESINS");

gs_file  <- paste(mark_prefix,".txt",sep="");
gene_id  <- paste(mark_prefix,"_gene_id",sep="");
gs_id    <- paste(mark_prefix,"_id",sep="");
gs_mod   <- paste(mark_prefix,"_id_gs_NW",sep="");
gs_Rdata <- paste(mark_prefix,".RData",sep="");

##
library(igraph)
library("RColorBrewer");
mypalette<-c(brewer.pal(8,"PuRd")[5],brewer.pal(8,"Greens")[5],brewer.pal(8,"Oranges")[4],brewer.pal(8,"Blues")[7],
		brewer.pal(6,"Greys")[4],brewer.pal(8,"Greens")[8],brewer.pal(8,"Purples")[7]);
##
if(!file.exists(gs_Rdata)){
	net2list <- function(x){
		x1 <- unique(x[,1]);
		x2 <- lapply(x1,function(y){x[which(x[,1] %in% y),2]})  
			names(x2) <- x1;x2
	}
##
	net_data <- unique(read.delim(net_file,stringsAsFactors=F,header=F));
	net_data <- net_data[which(net_data[,1] != net_data[,2]),]
		g1   <- graph.data.frame(net_data,directed=F) ## PPI
		net_g <- V(g1)$name;
	print("Finish read in network files")
		gs_data  <- unique(read.delim(gs_file,stringsAsFactors=F,header=F)); 
	gs_list <- net2list(gs_data[,2:1]); gs_list <- gs_list[unlist(lapply(gs_list,length)>=5)]
	names(gs_list) <- gsub("(GO........).*","\\1",names(gs_list));
		print("Finish read in gs files")
# 
		gene_id <- read.table(gene_id,stringsAsFactors=F)[,1]
		gs_id <- read.delim(gs_id,stringsAsFactors=F,header=F)[,1]
		gs_mod_data <- as.matrix(read.table(file=gs_mod)); ## three columns
		print(summary(as.numeric(gs_mod_data[,3])));
	rownames(gs_mod_data) <- paste(gene_id[gs_mod_data[,1]],gsub("(GO........).*","\\1",gs_id[gs_mod_data[,2]]),sep="~");
	gs_mod_data[,3] <- log10(2+gs_mod_data[,3]); 
	print(summary(as.numeric(gs_mod_data[,3])));
## 
	print("Finish read in bg files")
		save(list=ls(),file=gs_Rdata);
}

##
load(gs_Rdata);
##########
source("src/plot_LEGO_func.R");
input_enrich <- read.delim(input_file,header=T,stringsAsFactors=F)[,c(1,4,8)] ## read in LEGO output
input_enrich[,2] <- signif(input_enrich[,2],3) ## p-value
int_g <- read.delim(int_file,stringsAsFactors=F,header=F)[,1];
##
output_file_pdf <- paste(output_file,".pdf",sep="");
pdf(output_file_pdf,width=30,height=30);par(mar=c(8,5,8,5));par(cex.main=5);
all_rt <- NULL;
for(gs in gs_for_plot){
	x  <- input_enrich[which(gsub("(GO........).*","\\1",input_enrich[,1]) == gs | input_enrich[,1] == gs),];pv <- x[2];
	r1 <- plot_LEGO(gs, pv, int_g, choose = "full", inter = "in", ii = "", net_g = net_g, gs_list = gs_list, net_data = net_data, gs_mod_data, min_size=0,local=0);
	print(r1);
	plot(r1,main=paste(gs,pv));
	rt1 <- gs_mod_data[intersect(paste(V(r1)$name,gs,sep="~"),rownames(gs_mod_data)),3];
	names(rt1) <- gsub("(.*)~.*","\\1",names(rt1)); rt1 <- sort(rt1,decreasing=TRUE); rt1 <- signif(rt1,4); 
	rt2 <- rep(".",length(rt1)); rt2[which(names(rt1) %in% int_g)] <- "Y"; 
	rt3 <- rep(".",length(rt1)); rt3[which(names(rt1) %in% gs_list[[gs]])] <- gs; 
	rt <- data.frame(Gene_name=names(rt1),Score=rt1,Gene_set=rt3,Interesting_gene_list_index=rt2);
	all_rt <- rbind(all_rt,rt);
}
dev.off();
output_file_tab <- paste(output_file,".txt",sep="");
write.table(all_rt,file=output_file_tab,sep="\t",col.names=T,row.names=F,quote=F);
##
output_file_pdf <- paste(output_file,"_local.pdf",sep="")
pdf(output_file_pdf,width=20,height=20);par(mar=c(5,5,5,5));par(cex.main=3)
for(gs in gs_for_plot){
	x  <- input_enrich[which(gsub("(GO........).*","\\1",input_enrich[,1]) == gs | input_enrich[,1] == gs),];pv <- x[2];
	r2 <- plot_LEGO(gs, pv, int_g, choose = "full", inter = "in", ii = "", net_g = net_g, gs_list = gs_list, net_data = net_data, gs_mod_data, min_size=0,local=1);
	E(r2)$width <- E(r2)$width*0.5;
	V(r2)$size <- 1.5^V(r2)$size*1.5;
	print(r2);
	plot(r2,main=paste(gs,pv),layout=layout.kamada.kawai);
}
dev.off()
## 



