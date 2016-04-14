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
#Gene_set	Z_score(n-LEGO)	P_value(n-LEGO)	Z_score(e-LEGO)	P_value(e-LEGO)	gsSize	gsSize_net	intSize	intSize_net	overlapSize	InterestingGeneListID
#2	-0.151035	0.560028	0.220164	0.412872	33	25	170	151	3	1

cal_func <- function(f1){
	#p_Edge <- f1[,5];	
	p_Edge <- pnorm(f1[,4],lower.tail=F)
	padj_Edge <- as.numeric(p.adjust(p_Edge,method=met))
	p_Edge <- signif(p_Edge,5); padj_Edge <- signif(padj_Edge,5)

	f_Edge <- cbind(f1[,c(1,4)],p_Edge,padj_Edge,f1[,6:10])
	r_Edge <- order(f_Edge[,2],decreasing=T)
	f_Edge <- f_Edge[r_Edge,]
	colnames(f_Edge) <- c("GeneSet_des","Z_score","Ori_p","Adj_p","gsSize","gsSize_net","intSize","intSize_net","OverlapSize")

	return(list(f_Edge=f_Edge))
}	
################
f0  <- read.table(input,header=T)
len <- length(f0[1,])
if(multi==1){
	all_Edge <- NULL;
	all_int <- unique(f0[,len])
	for(i in all_int){
		uu <- which(f0[,len]==i); f1 <- f0[uu,]
		re <- cal_func(f1);	
		all_Edge <- rbind(all_Edge,cbind(re$f_Edge,InterestingGeneListID=rep(i,length(uu))))
	}
}else{
	re <- cal_func(f0);	
	all_Edge <- re$f_Edge
}
all_Edge <- all_Edge[which(all_Edge[,4]<= thre),]
write.table(all_Edge,file=paste(tmp,"_EdgeResults.txt",sep=""),quote=F,sep="\t",row.names=F)

