##
#Gene_set    Z_score(n-LEGO) Z_score(e-LEGO) gsSize  intSize overlapSize
#1   0.032200    0.111059    136 67  64

Args <- commandArgs(T)
input   <- Args[1]
thre <- as.numeric(Args[2]) ## p-value cutoff
met  <- Args[3] ## adjusted method
tmp  <- Args[4] ## output file
multi <- as.numeric(Args[5]) ## 
########################### functions
##
cal_func <- function(f1){
	p_Edge <- pnorm(f1[,3],lower.tail=F)
	padj_Edge <- as.numeric(p.adjust(p_Edge,method=met))
	p_Edge <- signif(p_Edge,5); padj_Edge <- signif(padj_Edge,5)

	p_Node <- pnorm(f1[,2],lower.tail=F)
	padj_Node <- as.numeric(p.adjust(p_Node,method=met))
	p_Node <- signif(p_Node,5); padj_Node <- signif(padj_Node,5)

	f_Edge <- cbind(f1[,c(1,3)],p_Edge,padj_Edge,f1[,4:8])
	r_Edge <- order(f_Edge[,2],decreasing=T)
	f_Edge <- f_Edge[r_Edge,]
	colnames(f_Edge) <- c("GeneSet_des","Z_score","Ori_p","Adj_p","gsSize","gsSize_net","intSize","intSize_net","OverlapSize")

	f_Node <- cbind(f1[,c(1,2)],p_Node,padj_Node,f1[,4:8])
	r_Node <- order(f_Node[,2],decreasing=T)
	f_Node <- f_Node[r_Node,]
	colnames(f_Node) <- c("GeneSet_des","Z_score","Ori_p","Adj_p","gsSize","gsSize_net","intSize","intSize_net","OverlapSize")

	return(list(f_Edge=f_Edge,f_Node=f_Node))
}	
################
f0  <- read.table(input,header=T)
len <- length(f0[1,])
if(multi==1){
	all_Edge <- NULL; all_Node <- NULL;
	all_int <- unique(f0[,len])
	for(i in all_int){
		uu <- which(f0[,len]==i); f1 <- f0[uu,]
		re <- cal_func(f1);	
		all_Edge <- rbind(all_Edge,cbind(re$f_Edge,InterestingGeneListID=rep(i,length(uu))))
		all_Node <- rbind(all_Node,cbind(re$f_Node,InterestingGeneListID=rep(i,length(uu))))
	}
}else{
	re <- cal_func(f0);	
	all_Edge <- re$f_Edge
	all_Node <- re$f_Node
}
all_Edge <- all_Edge[which(all_Edge[,4]< thre),]
all_Node <- all_Node[which(all_Node[,4]< thre),]
write.table(all_Edge,file=paste(tmp,"_EdgeResults.txt",sep=""),quote=F,sep="\t",row.names=F)
write.table(all_Node,file=paste(tmp,"_NodeResults.txt",sep=""),quote=F,sep="\t",row.names=F)

