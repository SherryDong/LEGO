#!~/bin/Rscript

LEGO_DIR = system('pwd',intern=TRUE)
Args <- commandArgs(T)
if(length(Args)<5){
	print("Usage: plot_LEGO.R <input interesting gene list file> <enriched output file> <output file prefix> <network file name> <network+gene set prefix,e.g: demo/GeneSet_human.txt_FC2_human> <CS or NW> <plot gene set name, can be comma saperated> \n");
	quit("no");
}
# system("R CMD SHLIB src/LEGO_draw.c")
dyn.load(paste(LEGO_DIR , '/src/LEGO_draw.so',sep=''))
library("popgraph")
##
int_file    <- Args[1]
input_file  <- Args[2] ## 
output_file <- Args[3]
net_file    <- Args[4]
mark_prefix <- Args[5]
cs_or_nw    <- Args[6];
gs_for_plot <- Args[7]
gs_for_plot <- unlist(strsplit(gs_for_plot,",")); 
print(gs_for_plot)

###for example:
# Rscript src/plot_LEGO.R  demo_server/autism_2.txt autism_2_out/autism_2.txt_id_MSigDB_1.0.out_LEGO.txt autism_2_out/autism_2_GABA /home/sherry/project/network/NETGO/LEGO_server_v2.0/data/FunCoup/FC3.0_H.sapiens_compact /home/sherry/project/network/NETGO/LEGO_server_v2.0/data/GeneSet/MSigDB_1.0/goa_human.res_FC3.0_H.sapiens_compact CS BIOCARTA_GABA_PATHWAY 

#int_file <- "demo/input.txt";
#input_file <- "demo/input.txt_id.out_e-LEGO.txt";
#output_file <- "demo/test_plot";
#net_file <- "demo/FC2_human.txt";
#mark_prefix <- "demo/GeneSet_human.txt_FC2_human";
#gs_for_plot <- c("REACTOME_G1_S_SPECIFIC_TRANSCRIPTION","REACTOME_KINESINS");

gs_file  <- mark_prefix;
gene_id  <- paste(mark_prefix,"_gene_id",sep="");
gs_id    <- paste(mark_prefix,"_id",sep="");
int_id   <- paste(int_file,"_id",sep="");
gs_id_gs <- paste(mark_prefix,"_id_gs",sep="");
net_id_file <- gene_id;
#net_id_file <- paste(net_file,"_id",sep="");
net_id_net_file <- paste(net_file,"_id_net",sep="");

if(cs_or_nw=="CS"){
	gs_mod   <- paste(mark_prefix,"_id_gs_CS",sep="");
}else{
	gs_mod   <- paste(mark_prefix,"_id_gs_NW",sep="");
}
##
library(igraph)
library("RColorBrewer");
mypalette<-c(brewer.pal(8,"PuRd")[5],brewer.pal(8,"Greens")[5],brewer.pal(8,"Oranges")[4],brewer.pal(8,"Blues")[7],
		brewer.pal(6,"Greys")[4],brewer.pal(8,"Greens")[8],brewer.pal(8,"Purples")[7]);
##
##########
source(paste(LEGO_DIR ,"/src/plot_LEGO_func.R", sep=''));

input_enrich <- read.delim(input_file,header=T,stringsAsFactors=F)[,c(1,2)]; ## read in LEGO output
input_enrich[,2] <- signif(input_enrich[,2],3); ## p-value
int_g_id <- read.delim(int_id,stringsAsFactors=F,header=F)[,1]; ## input gene list id
## get gs id
gs2id_tmp <- read.delim(gs_id,stringsAsFactors=F,header=F); gs2id <- gs2id_tmp[,2]; names(gs2id) <- gs2id_tmp[,1]; 
id2gene_tmp <- read.delim(net_id_file,stringsAsFactors=F,header=F); id2gene <- id2gene_tmp[,1]; names(id2gene) <- as.character(id2gene_tmp[,2]);
##
##
input_html = paste(LEGO_DIR,"/www/html/LEGO_template.html",sep=""); 
output_file_pdf <- paste(output_file,".pdf",sep="");
output_file_gml <- paste(output_file,".graphml",sep="");
output_file_html <- paste(output_file,".html",sep="");
pdf(output_file_pdf,width=30,height=30);par(mar=c(8,5,8,5));par(cex.main=5);
for(gs in gs_for_plot){
	print(gs);
	use_gs_id <- gs2id[gs]; ## get gene set id 

	out1 <- .Call("LEGO_draw_gene_set",c(gs_id_gs,use_gs_id),package="LEGO_draw"); ## extract gene id in the gene set 
	all_use_gene_id <- unique(c(out1,int_g_id)); ## all input and gene set genes 
	para <- c(net_id_file,net_id_net_file,length(all_use_gene_id),all_use_gene_id);
	out2 <- .Call("LEGO_draw_net_file",para,package="LEGO_draw"); ## extract sub network
	sub_net <- matrix(out2,ncol=2,byrow=TRUE);
	para <- c(net_id_file,gs_mod,use_gs_id,length(all_use_gene_id),all_use_gene_id); ## file, gene set id, gene id length, gene id	
	out3 <- .Call("LEGO_draw_gs_mod",para,package="LEGO_draw"); ## input gene set name + use gene id 	
	gs_mod_data <- matrix(out3,ncol=2,byrow=TRUE); rownames(gs_mod_data) <- gs_mod_data[,1]; 
	gs_mod_data[,2] <- log10(2+gs_mod_data[,2]); 
	x  <- input_enrich[which(gsub("(GO........).*","\\1",input_enrich[,1]) == gs | input_enrich[,1] == gs),];pv <- x[2];
	r1 <- plot_LEGO(gs, pv, net_g = sub_net, all_use_gene_id, int_g=setdiff(int_g_id,out1),gs_g=setdiff(out1,int_g_id),ov=intersect(out1,int_g_id),gs_mod_data,min_size=0,local=0,id2gene);

	print(r1);plot(r1,main=paste(gs,pv));
	tmp_file_name <- graph2xml(r1, input_html, output_file_html);
	if(gs == gs_for_plot[1]){
		write.graph(r1, file = output_file_gml, format = "graphml")
	}else{
		write.graph(r1, file = output_file_gml, format = "graphml",append=TRUE);
	}
}
dev.off();
##
#output_file_pdf <- paste(output_file,"_local.pdf",sep="")
#output_file_gml <- paste(output_file,"_local.gml",sep="")
# to_json(graph, file)
#pdf(output_file_pdf,width=20,height=20);par(mar=c(5,5,5,5));par(cex.main=3)
#for(gs in gs_for_plot){
#	print(gs);x  <- input_enrich[which(gsub("(GO........).*","\\1",input_enrich[,1]) == gs | input_enrich[,1] == gs),];pv <- x[2];
#	r2 <- plot_LEGO(gs, pv, net_g = sub_net, all_use_gene_id, int_g=setdiff(int_g_id,out1),gs_g=setdiff(out1,int_g_id),ov=intersect(out1,int_g_id),gs_mod_data,min_size=0,local=1,id2gene);
#	E(r2)$width <- E(r2)$width*0.5; V(r2)$size <- 1.5^V(r2)$size*1.5;
#	print(r2); plot(r2,main=paste(gs,pv),layout=layout.kamada.kawai);
#	if(gs == gs_for_plot[1]){
#		write.graph(r2, file = output_file_gml, format = "graphml")
#		e2 <- as.popgraph(r2); E(e2)$weight <- 1;
#	}else{
#		write.graph(r2, file = output_file_gml, format = "graphml",append=TRUE);
#		e2 <- as.popgraph(r2); E(e2)$weight <- 1;
#	}
#}
#dev.off()
## 

