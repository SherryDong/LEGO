##
#Gene_set    Z_score(n-LEGO) Z_score(e-LEGO) gsSize  intSize overlapSize
#1   0.032200    0.111059    136 67  64
library(multicore)
Args <- commandArgs(T)
input   <- Args[1]
thre <- as.numeric(Args[2]) ## p-value cutoff
met  <- Args[3] ## adjusted method
tmp  <- Args[4] ## output file
multi <- as.numeric(Args[5]) ## 
mid_file_path <- Args[6];
###
p_thre <- 0.05;
if(thre < p_thre){
	p_thre <- thre;
}
library(fExtremes);
F_gpd <- function(all_s){
	all_s <- sort(all_s,decreasing=TRUE);
	while(test_gpd(all_s)<0.05 & length(all_s)>100){
		#print(c(length(all_s),test_gpd(all_s)));
		all_s <- all_s[1:(length(all_s)-10)];
	}   
	fit_d <-  gpdFit(all_s,u=0);
	return(fit_d);
}
test_gpd <- function(all_s){
	fit_d <-  gpdFit(all_s,u=0);
	sim_d <-  as.numeric(gpdSim(model = list(beta=fit_d@fit$par.ests["beta"],mu=0,xi=fit_d@fit$par.ests["xi"]), 
				n = length(all_s)));
	pv<-ks.test(all_s,sim_d)$p.value;return(pv);
}    
pv_gpd <- function(data_rank,ob,p_thre){
	numT=250;
	data_rank <- sort(data_rank,decreasing=TRUE);
	ori_p <- (1+length(which(data_rank>=ob)))/length(data_rank);
	if(ori_p>1) ori_p <- 1;
	if(ori_p>p_thre){return(ori_p)};
	t_score <- (data_rank[numT]+data_rank[numT+1])/2;
	fit_d <- try(F_gpd(c(ob,data_rank[1:numT])-t_score));
	if(class(fit_d)=="try-error"){return(ori_p);}
	pv <- pgpd(ob-t_score,mu=0,beta=fit_d@fit$par.ests["beta"],xi=fit_d@fit$par.ests["xi"]);
	adj_p <- ori_p*(1-pv);
	return(adj_p);
}
########################### functions
##
#Gene_set	Z_score(n-LEGO)	P_value(n-LEGO)	Z_score(e-LEGO)	P_value(e-LEGO)	gsSize	gsSize_net	intSize	intSize_net	overlapSize	InterestingGeneListID
#2	-0.151035	0.560028	0.220164	0.412872	33	100	170	151	3	1
cal_func <- function(f1){
	p_Edge <- unlist(mclapply(c(1:length(f1[,1])),function(x){
				gs = f1[x,1];ob = f1[x,4];
				w1 <- which(all_gs == gs); 
				w2 <- c(0:(perm_times-1))*length(all_gs)+w1;
				data_rank <- input_info[w2]; 
				return(pv_gpd(data_rank,ob,p_thre));
			},mc.cores=15))	
	padj_Edge <- as.numeric(p.adjust(p_Edge,method=met))
	p_Edge <- signif(p_Edge,2); padj_Edge <- signif(padj_Edge,2)
	
	f_Edge <- cbind(f1[,1],p_Edge,padj_Edge,f1[,6:10])
	r_Edge <- order(f_Edge[,2],decreasing=FALSE)
	f_Edge <- f_Edge[r_Edge,]
	colnames(f_Edge) <- c("GeneSet_des","Ori_p","Adj_p","gsSize","gsSize_net","intSize","intSize_net","OverlapSize")

	return(list(f_Edge=f_Edge))
}	
################
f0  <- read.table(input,header=T)
len <- length(f0[1,])
if(multi==1){
	all_Edge <- NULL; 
	all_int <- unique(f0[,len])
	for(i in all_int){
		uu <- which(f0[,len]==i); f1 <- f0[uu,]; 
		intSize <- f1[1,"intSize_net"];
		file_perm <- paste(mid_file_path,intSize,".RData",sep="");print(file_perm);load(file_perm);	
		perm_times <- round(length(input_info)/length(all_gs));print(perm_times);
		re <- cal_func(f1);	
		all_Edge <- rbind(all_Edge,cbind(re$f_Edge,InterestingGeneListID=rep(i,length(uu))))
	}
}else{
	f1 <- f0;
	intSize <- f1[1,"intSize_net"];
	file_perm <- paste(mid_file_path,intSize,".RData",sep="");print(file_perm);load(file_perm);	
	perm_times <- round(length(input_info)/length(all_gs));print(perm_times);
	re <- cal_func(f1);	
	all_Edge <- re$f_Edge
}
all_Edge <- all_Edge[which(all_Edge[,3]<= thre),]
write.table(all_Edge,file=paste(tmp,"_EdgeResults.txt",sep=""),quote=F,sep="\t",row.names=F)

