/* default parameters for basic graph task */
#ifndef _io
#define _io

#include <stdio.h>
#include <stdlib.h>	//itoa	bsearch

#include <string.h> //strcpy ...

#include <math.h> //log ...
//-lm when compile

#include <stdarg.h>	//va_list
#include <errno.h> //errno, strerror
#include <unistd.h> //getopt
#include <time.h> //time

#include <sys/stat.h>	//mkdir
#include <sys/types.h>

#include <signal.h>	//fork
#include <wait.h>	//wait

//#include <limits.h>
#include <float.h>	//DBL_MIN
#include <limits.h> // INT_MAX

#define FALSE 0
#define TRUE 1
#define False 0
#define True 1
#define F 0
#define T 1

//array macros
#define MAX_NUM_SIZE 20
#define MAX_FILE_NAME_LENGTH 1000
#define SINGLE_PRAR_LENGTH 50
#define COMMAND_LINE_STRING_LENGTH 1000

#define MAX_LINE_SIZE 100000
#define MAX_NAME_SIZE 100
#define MAX_NODE_SIZE 100000
#define MAX_EDGE_SIZE 5000000
#define MAX_GS_SIZE 20000
#define MAX_IG_SIZE 100000
#define INF 0x0fffffff  

#define DELIMITER '\t'
//macro as functions
#define MAX(A, B) A<B?B:A
#define MIN(A, B) A<B?A:B
typedef _Bool boolean;
typedef double NUM;
typedef int COUNT;
// net structure
struct ArcNode // edge
{
	int adjvex; // index for another node in edge
	double weight; // edge weight
	struct ArcNode *nextarc; // point to next edge
};
struct VNode // node
{
	int data; // node index
	struct ArcNode *head; // point
};
struct LGraph // graph structure
{
	struct VNode vertexs[MAX_NODE_SIZE]; // array for node
	int vexnum,arcnum; // number of node, number of edge
};
struct GS_INFO
{
	int size; // size
	int size_net; // size
	int ov;
	double n_mean;
	double n_sd;
	double e_mean;
	double e_sd;
	double n_ES;
	double e_ES;
};
struct GS_INFO_BG
{
	int size; // size
	int size_net; // size
	int ov;
	double n_mean;
	double n_sd;
	double n_ES;
	double n_ES_bg;
	
	double e_mean;
	double e_sd;
	double e_ES;
	double e_ES_bg;
};
#define PI 3.14159265359
double z2p(double x){
	double b0 = 0.2316419;
	double b1 = 0.319381530;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double t,pdf,cdf,p;
	t = 1/(1+b0*x);
	pdf = exp(-x*x/2)/sqrt(2*PI);
	cdf = 1-pdf*(b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5));
	p = 1-cdf;
	return(p);
}
int wc_l(char *file){
	FILE *fp;
	char line[MAX_LINE_SIZE];
	COUNT line_num;
	if ((fp = fopen(file, "r")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s;\n",
			file, strerror(errno) );
		return -1;
	}
	line_num=0;
	while (fgets(line, sizeof(line), fp) != NULL 
		&& line[0] != '\r' && line[0] != '\0') {
		line_num++;
	}
	fclose(fp);
	return line_num;
}
// get gs info
int GetGS_info(char *input_gs_file,int *nodeNum,int *gsNum){
	FILE *fp;
	int i,j,d1,d2,tmp1,tmp2;
	if ((fp = fopen(input_gs_file, "r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_gs_file, strerror(errno) );
		return -1;
	}
	tmp1=0;tmp2=0;
	while(fscanf(fp,"%d\t%d",&d1,&d2)!=EOF){
		tmp1=(tmp1>d1)?tmp1:d1;
		tmp2=(tmp2>d2)?tmp2:d2;
	}
	(*nodeNum) = tmp1;
	(*gsNum) = tmp2;
	printf("In total %d genes and %d gene sets\n",tmp1,tmp2);
	fclose(fp);
}
// read GS from file
int ReadGS_file(char *input_gs_file,int nodeNum,int vexnum,boolean ***node_gs_ptr,int **gsSize,int **gsSize_net,int gsNum){
	FILE *fp;
	int i,j,d1,d2,tmp;
	if ((fp = fopen(input_gs_file, "r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_gs_file, strerror(errno) );
		return -1;
	}
	for(i=1;i<=gsNum;i++){
		(*gsSize)[i]=0;
		(*gsSize_net)[i]=0;
		for(j=1;j<=nodeNum;j++){
			(*node_gs_ptr)[j][i] = 0;
		}
	}
	fp = fopen(input_gs_file, "r");
	while(fscanf(fp,"%d\t%d",&d1,&d2)!=EOF){
		(*node_gs_ptr)[d1][d2] = 1;
	}
	fclose(fp);
	printf("Successful read in %s, total gene set number: %d\n",input_gs_file,gsNum);
	// get gene set size in network
	for(i=1;i<=gsNum;i++){
		for(j=1;j<=nodeNum;j++){
			if((*node_gs_ptr)[j][i] ==1){
				(*gsSize)[i] ++;	
				if(j<=vexnum){
					(*gsSize_net)[i] ++;
				}
			}
		}
		//printf("Size for gene set in network %d: %d\n",i,(*gsSize_net)[i]);
	}
	return(1);
}
// read multiple IG from file
int ReadIG_mul_file(char *input_int_file,char *input_GS_file,boolean ***intG,int **intGNum,int *intNum,int **intGNum_net){	
	FILE *fp1,*fp2;
	int d1,id,i,j,m1,nodeNum,vexnum,m4,m5,m6;
	if ((fp1 = fopen(input_int_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_int_file, strerror(errno));
		return -1;
	}
	if ((fp2 = fopen(input_GS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_GS_file, strerror(errno));
		return -1;
	}
	while(fscanf(fp2,"%d\t%d\t%d\t%d\t%d\t%d\n",&m1,&nodeNum,&vexnum,&m4,&m5,&m6)!=EOF){
		break;	
	}
	for(i=0;i<=MAX_IG_SIZE;i++){
		(*intG)[i]  = (boolean *)malloc((nodeNum+1)*sizeof(boolean));
	}   
	(*intNum)=0;
	while(fscanf(fp1,"%d\t%d\n",&d1,&id)!=EOF){
		(*intNum)=((*intNum)>id)?(*intNum):id;
	}
	fclose(fp1);
	for(i=1;i<=(*intNum);i++){
		(*intGNum)[i]=0;
		(*intGNum_net)[i]=0;
		for(j=1;j<=vexnum;j++){
			(*intG)[i][j]=0;
		}
	}
	fp1 = fopen(input_int_file,"r");
	while(fscanf(fp1,"%d\t%d\n",&d1,&id)!=EOF){
		(*intG)[id][d1]=1;(*intGNum)[id]++;
		if(d1<=vexnum){(*intGNum_net)[id]++;}
	}
	fclose(fp1);
	printf("Finish read in interesting gene list !\n");
	return(1);
}
// Read IG from file
int ReadIG_int_file (char *input_int_file,char *input_GS_file,boolean **intG, int *intSize, int *intSize_net){	
	FILE *fp1,*fp2;
	int d1,m1,m2,vexnum,m4,m5,m6;
	if ((fp1 = fopen(input_int_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_int_file, strerror(errno));
		return -1;
	}
	if ((fp2 = fopen(input_GS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_GS_file, strerror(errno));
		return -1;
	}
	while(fscanf(fp2,"%d\t%d\t%d\t%d\t%d\t%d\n",&m1,&m2,&vexnum,&m4,&m5,&m6)!=EOF){
		break;	
	}
	(*intSize)=0; (*intSize_net)=0;
	while(fscanf(fp1,"%d\n",&d1)!=EOF){
		(*intG)[d1]=1;
		(*intSize)++;
		if(d1<=vexnum){(*intSize_net)++;}
	}
	printf("Finish read in interesting gene list !\n");
	return(1);
}
// output_NW_file,output_GS_file,output_CS_file,lg,node_wei,gsNum,node_gs,gsSize
int Out_NW_CS_file(char *output_NW_file,char *output_GS_file, char *output_CS_file, struct LGraph lg,double **node_wei,int gsNum,boolean **node_gs,int *gsSize, int *gsSize_net,int nodeNum){
	FILE *fp1, *fp2, *fp3; 
	double wei,wei_e,*edgeS, bg_sum, bg_v1, avg, sd;
	int i,j,d1,d2;
	struct ArcNode *pi; // edge pointer
	pi=(struct ArcNode *)malloc(sizeof(struct ArcNode));
	edgeS=(double*)malloc((nodeNum+1)*sizeof(double));
	if ((fp1 = fopen(output_NW_file,"w")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",output_NW_file, strerror(errno));
		return -1;
	}
	if ((fp2 = fopen(output_CS_file,"w")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",output_GS_file, strerror(errno));
		return -1;
	}
	if ((fp3 = fopen(output_GS_file,"w")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",output_CS_file, strerror(errno));
		return -1;
	}
	fprintf(fp3,"0\t%d\t%d\t%d\t0\t0\n",nodeNum,lg.vexnum,lg.arcnum);
	for(i=1;i<=gsNum;i++){ // for each gene set
		// n-LEGO
		for(j=1;j<=nodeNum;j++){ // for each gene in <geneset>_gene_id
			if(node_gs[j][i]==1){ // belong to gs
				if(j <= lg.vexnum){
					edgeS[j]=node_wei[j][i]; // belong to network
				}else{
					edgeS[j]=1; // does not belong to network
				}
			}else{
				edgeS[j]=0;	
			}
		}
		bg_sum=0; bg_v1=0;
		for(j=1;j<=nodeNum;j++){ // for each gene, print out NW
			if(edgeS[j]>0){
				fprintf(fp1,"%d\t%d\t%f\n",j,i,edgeS[j]);
			}
			bg_sum+=edgeS[j];bg_v1+=edgeS[j]*edgeS[j];
			edgeS[j]=0; // default = 0 for e-LEGO
		}
		avg = bg_sum/nodeNum;
		sd = sqrt(bg_v1/nodeNum-avg*avg); 
		fprintf(fp3,"1\t%d\t%d\t%d\t%f\t%f\n",i,gsSize[i],gsSize_net[i],avg,sd);
		// e-LEGO
		for(j=1;j<=lg.vexnum;j++){ // for each gene in network
			if(node_gs[j][i]==1){ // belong to gs
				wei_e=node_wei[j][i]*node_wei[j][i]; // self weight
				edgeS[j]+=wei_e;
				pi=lg.vertexs[j-1].head; // for each gene belong to gs
				d1=j;
				while(pi!=NULL){
					wei=pi->weight;
					d2=(pi->adjvex)+1; // interaction ID
					wei_e=wei*node_wei[d1][i]*node_wei[d2][i];
					edgeS[d2]+=wei_e;
					pi=pi->nextarc;
				}
			}
		}
		bg_sum=0; bg_v1=0;
		for(j=1;j<=lg.vexnum;j++){ // for each gene
			if(edgeS[j]>0){
				fprintf(fp2,"%d\t%d\t%f\n",j,i,edgeS[j]);
			}
			bg_sum+=edgeS[j];bg_v1+=edgeS[j]*edgeS[j];
		}
		avg = bg_sum/lg.vexnum;
		sd = sqrt(bg_v1/lg.vexnum-avg*avg); 
		fprintf(fp3,"2\t%d\t%d\t%d\t%f\t%f\n",i,gsSize[i],gsSize_net[i],avg,sd);
	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	return 1;
}
// output results
LEGOout_int_file(boolean *intG,int intSize,int intSize_net,int gsNum,int geneset_min, int geneset_max,char *output_ES_file,char *input_NW_file,char *input_CS_file,char *input_GS_file){
	FILE *fp1,*fp2,*fp3,*fp4;
	int i,mm,gs_ID,gene_ID,gs_size,gs_size_net,bgSize;
	double avg,sd,wei,Z_score_n,Z_score_e,P_score_n,P_score_e;
	boolean *use_gs;
	struct GS_INFO *gs;
	gs   = (struct GS_INFO *)malloc((gsNum+1)*sizeof(struct GS_INFO));
	use_gs = (boolean *)malloc((gsNum+1)*sizeof(boolean));
	for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
		use_gs[gs_ID]=0;
	}
	if ((fp1 = fopen(output_ES_file,"w")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",output_ES_file, strerror(errno));
		return -1;
	}
	if ((fp2 = fopen(input_NW_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_NW_file, strerror(errno));
		return -1;
	}
	if ((fp3 = fopen(input_CS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_CS_file, strerror(errno));
		return -1;
	}
	if ((fp4 = fopen(input_GS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_GS_file, strerror(errno));
		return -1;
	}
	// read in GS info file 
	// 0	18716	17083	4353383	0	0
	// 1	1	16	14	0.063120	3.386597
	// e.g 1       1       52      0.137187        5.216251	
	// 1   1015    10  4   0.003192    0.223718
	i = 0;
	while(fscanf(fp4,"%d\t%d\t%d\t%d\t%lf\t%lf\n",&mm,&gs_ID,&gs_size,&gs_size_net,&avg,&sd)!=EOF){
		i++;
		if(i==1)
			bgSize = gs_size;
		if(gs_size < geneset_min || gs_size > geneset_max){ // size selection
			continue;
		}
		if(mm==1){ // n-LEGO
			use_gs[gs_ID] = 1; // use gs
			gs[gs_ID].size = gs_size;
			gs[gs_ID].size_net = gs_size_net;
			gs[gs_ID].n_mean = avg;
			gs[gs_ID].n_sd = sd;
			gs[gs_ID].n_ES = 0;
			gs[gs_ID].ov   = 0;
		}else if(mm==2){
			gs[gs_ID].e_mean = avg;
			gs[gs_ID].e_sd = sd;
			gs[gs_ID].e_ES = 0;
		} 
	}		
	fclose(fp4);
	printf("finish read in %s\n",input_GS_file);
	// read in NW 
	// 786	1	83.742686
	while(fscanf(fp2,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			if(intG[gene_ID]==1){ // belong to int
				gs[gs_ID].n_ES += wei;						
				gs[gs_ID].ov ++; // intersect all
			}	
		}
	}
	fclose(fp2);
	printf("finish read in %s\n",input_NW_file);
	// read in CS
	while(fscanf(fp3,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			if(intG[gene_ID]==1){ // belong to int
				gs[gs_ID].e_ES += wei;						
			}	
		}
	}
	fclose(fp3);
	printf("finish read in %s\n",input_CS_file);
	// output
	fprintf(fp1,"Gene_set\tZ_score(n-LEGO)\tP_value(n-LEGO)\tZ_score(e-LEGO)\tP_value(e-LEGO)\tgsSize\tgsSize_net\tintSize\tintSize_net\toverlapSize\n");
	double scale_val;
	for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
		if(use_gs[gs_ID]){ // use gs
			gs[gs_ID].n_ES = gs[gs_ID].n_ES/intSize;	
			scale_val = 1/sqrt(1-(double)intSize/(double)bgSize);
			Z_score_n = scale_val*(gs[gs_ID].n_ES-gs[gs_ID].n_mean)/(gs[gs_ID].n_sd/sqrt(intSize));
			P_score_n = z2p(Z_score_n);
			gs[gs_ID].e_ES = gs[gs_ID].e_ES/intSize_net;	
			scale_val = 1/sqrt(1-(double)intSize_net/(double)bgSize);
			Z_score_e = scale_val*(gs[gs_ID].e_ES-gs[gs_ID].e_mean)/(gs[gs_ID].e_sd/sqrt(intSize_net));
			P_score_e = z2p(Z_score_e);
			if(P_score_n>1 || P_score_n <0)
				P_score_n = 1;
			if(P_score_e>1 || P_score_e <0)
				P_score_e = 1;
			fprintf(fp1,"%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\n",gs_ID,Z_score_n,P_score_n,Z_score_e,P_score_e,gs[gs_ID].size,gs[gs_ID].size_net,intSize,intSize_net,gs[gs_ID].ov);
		}
	}
	fclose(fp1);
	return 1;
}
// output multi results
// LEGOout_int_mul_file((intG,intGNum,intNum,gsNum,geneset_min,geneset_max,output_ES_file,input_NW_file,input_CS_file,input_GS_file);)
LEGOout_int_mul_file(boolean **intG,int *intGNum,int *intGNum_net, int intNum, int gsNum,int geneset_min, int geneset_max,char *output_ES_file,char *input_NW_file,char *input_CS_file,char *input_GS_file){
	FILE *fp1,*fp2,*fp3,*fp4;
	int i,mm,gs_ID,gene_ID,gs_size,gs_size_net;
	double avg,sd,wei,Z_score_n,Z_score_e,P_score_n,P_score_e;
	boolean *use_gs;
	struct GS_INFO **gs;
	gs  = (struct GS_INFO **)malloc((intNum+1)*sizeof(struct GS_INFO *));
	for(i=1;i<=intNum;i++){
		gs[i]   = (struct GS_INFO *)malloc((gsNum+1)*sizeof(struct GS_INFO));
	}
	use_gs = (boolean *)malloc((gsNum+1)*sizeof(boolean));
	for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
		use_gs[gs_ID]=0;
	}
	if ((fp1 = fopen(output_ES_file,"w")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",output_ES_file, strerror(errno));
		return -1;
	}
	if ((fp2 = fopen(input_NW_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_NW_file, strerror(errno));
		return -1;
	}
	if ((fp3 = fopen(input_CS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_CS_file, strerror(errno));
		return -1;
	}
	if ((fp4 = fopen(input_GS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_GS_file, strerror(errno));
		return -1;
	}
	// read in GS info file 
	// e.g 1       1       52      0.137187        5.216251	
	int bgSize;
	i = 0;
	while(fscanf(fp4,"%d\t%d\t%d\t%d\t%lf\t%lf\n",&mm,&gs_ID,&gs_size,&gs_size_net,&avg,&sd)!=EOF){
		i++;
		if(i==1)
			bgSize = gs_size;
		if(gs_size < geneset_min || gs_size > geneset_max){ // size selection
			continue;
		}
		for(i=1;i<=intNum;i++){
			if(mm==1){ // n-LEGO
				use_gs[gs_ID] = 1; // use gs
				gs[i][gs_ID].size = gs_size;
				gs[i][gs_ID].size_net = gs_size_net;
				gs[i][gs_ID].n_mean = avg;
				gs[i][gs_ID].n_sd = sd;
				gs[i][gs_ID].n_ES = 0;
				gs[i][gs_ID].ov = 0;
			}
			if(mm==2){
				gs[i][gs_ID].e_mean = avg;
				gs[i][gs_ID].e_sd = sd;
				gs[i][gs_ID].e_ES = 0;
			}
		}
	}		
	fclose(fp4);
	printf("finish read in %s\n",input_GS_file);
	// read in NW 
	// 786	1	83.742686
	while(fscanf(fp2,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			for(i=1;i<=intNum;i++){
				if(intG[i][gene_ID]==1){ // belong to int
					gs[i][gs_ID].n_ES += wei;						
					gs[i][gs_ID].ov ++;
				}
			}
		}
	}
	fclose(fp2);
	printf("finish read in %s\n",input_NW_file);
	// read in CS
	while(fscanf(fp3,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			for(i=1;i<=intNum;i++){
				if(intG[i][gene_ID]==1){ // belong to int
					gs[i][gs_ID].e_ES += wei;						
				}
			}	
		}
	}
	fclose(fp3);
	printf("finish read in %s\n",input_CS_file);
	// output
	double scale_val;
	fprintf(fp1,"Gene_set\tZ_score(n-LEGO)\tP_value(n-LEGO)\tZ_score(e-LEGO)\tP_value(e-LEGO)\tgsSize\tgsSize_net\tintSize\tintSize_net\toverlapSize\tInterestingGeneListID\n");
	for(i=1;i<=intNum;i++){
		for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
			if(use_gs[gs_ID]){ // use gs
	//			printf("%d\t%d\n",gs_ID,use_gs[gs_ID]);
				gs[i][gs_ID].n_ES = gs[i][gs_ID].n_ES/intGNum[i];	
				scale_val = 1/sqrt(1-(double)intGNum[i]/(double)bgSize);
				Z_score_n = scale_val*(gs[i][gs_ID].n_ES-gs[i][gs_ID].n_mean)/(gs[i][gs_ID].n_sd/sqrt(intGNum[i]));
				P_score_n = z2p(Z_score_n);
				gs[i][gs_ID].e_ES = gs[i][gs_ID].e_ES/intGNum_net[i];	
				scale_val = 1/sqrt(1-(double)intGNum_net[i]/(double)bgSize);
				Z_score_e = scale_val*(gs[i][gs_ID].e_ES-gs[i][gs_ID].e_mean)/(gs[i][gs_ID].e_sd/sqrt(intGNum_net[i]));
				P_score_e = z2p(Z_score_e);
				if(P_score_n>1 || P_score_n < 0)
					P_score_n = 1;
				if(P_score_e>1 || P_score_e < 0)
					P_score_e = 1;
				fprintf(fp1,"%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\n",gs_ID,Z_score_n,P_score_n,Z_score_e,P_score_e,gs[i][gs_ID].size,gs[i][gs_ID].size_net,intGNum[i],intGNum_net[i],gs[i][gs_ID].ov,i);
			}
		}
	}
	fclose(fp1);
	return 1;
}
// output with bg file
// output results
LEGOout_int_bg_file(boolean *intG,boolean *bgG,int intSize,int intSize_net,int bgSize, int bgSize_net, int gsNum,int geneset_min, int geneset_max,char *output_ES_file,char *input_NW_file,char *input_CS_file,char *input_GS_file){
	FILE *fp1,*fp2,*fp3,*fp4;
	int mm,gs_ID,gene_ID,gs_size,gs_size_net;
	double avg,sd,wei,Z_score_n,Z_score_e,P_score_n,P_score_e;
	boolean *use_gs;
	struct GS_INFO_BG *gs;
	gs   = (struct GS_INFO_BG *)malloc((gsNum+1)*sizeof(struct GS_INFO_BG));
	use_gs = (boolean *)malloc((gsNum+1)*sizeof(boolean));
	for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
		use_gs[gs_ID]=0;
	}
	if ((fp1 = fopen(output_ES_file,"w")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",output_ES_file, strerror(errno));
		return -1;
	}
	if ((fp2 = fopen(input_NW_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_NW_file, strerror(errno));
		return -1;
	}
	if ((fp3 = fopen(input_CS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_CS_file, strerror(errno));
		return -1;
	}
	if ((fp4 = fopen(input_GS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_GS_file, strerror(errno));
		return -1;
	}
	// read in GS info file 
	// e.g 1       1       52      0.137187        5.216251	
	while(fscanf(fp4,"%d\t%d\t%d\t%d\t%lf\t%lf\n",&mm,&gs_ID,&gs_size,&gs_size_net,&avg,&sd)!=EOF){
		if(gs_size < geneset_min || gs_size > geneset_max){ // size selection
			continue;
		}
		if(mm==1){ // n-LEGO
			use_gs[gs_ID] = 1; // use gs
			gs[gs_ID].size = gs_size;
			gs[gs_ID].size_net = gs_size_net;
			gs[gs_ID].n_mean = avg;
			gs[gs_ID].n_sd = sd;
			gs[gs_ID].n_ES = 0;
			gs[gs_ID].n_ES_bg = 0;
			gs[gs_ID].ov   = 0;
		}else if(mm==2){
			gs[gs_ID].e_mean = avg;
			gs[gs_ID].e_sd = sd;
			gs[gs_ID].e_ES = 0;
			gs[gs_ID].e_ES_bg = 0;
		} 
	}		
	fclose(fp4);
	printf("finish read in %s\n",input_GS_file);
	// read in NW 
	// 786	1	83.742686
	while(fscanf(fp2,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			if(intG[gene_ID]==1){ // belong to int
				gs[gs_ID].n_ES += wei;						
				gs[gs_ID].ov ++; // intersect all
			}	
			if(bgG[gene_ID]==1){ // belong to int
				gs[gs_ID].n_ES_bg += wei;	
			}
		}
	}
	fclose(fp2);
	printf("finish read in %s\n",input_NW_file);
	// read in CS
	while(fscanf(fp3,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			if(intG[gene_ID]==1){ // belong to int
				gs[gs_ID].e_ES += wei;						
			}	
			if(bgG[gene_ID]==1){ // belong to int
				gs[gs_ID].e_ES_bg += wei;	
			}
		}
	}
	fclose(fp3);
	printf("finish read in %s\n",input_CS_file);
	// output
	double scale_val;
	fprintf(fp1,"Gene_set\tZ_score(n-LEGO)\tP_value(n-LEGO)\tZ_score(e-LEGO)\tP_value(e-LEGO)\tgsSize\tgsSize_net\tintSize\tintSize_net\toverlapSize\n");
	for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
		if(use_gs[gs_ID]){ // use gs
			gs[gs_ID].n_ES = gs[gs_ID].n_ES/intSize;	
			gs[gs_ID].n_ES_bg = gs[gs_ID].n_ES_bg/bgSize;	
			scale_val = 1/sqrt(1-(double)intSize/(double)bgSize);
			Z_score_n = scale_val*(gs[gs_ID].n_ES-gs[gs_ID].n_ES_bg)/(gs[gs_ID].n_sd*sqrt(1/(double)intSize+1/(double)bgSize));
			P_score_n = z2p(Z_score_n);
			gs[gs_ID].e_ES = gs[gs_ID].e_ES/intSize_net;	
			gs[gs_ID].e_ES_bg = gs[gs_ID].e_ES_bg/bgSize_net;	
			scale_val = 1/sqrt(1-(double)intSize_net/(double)bgSize);
			Z_score_e = scale_val*(gs[gs_ID].e_ES-gs[gs_ID].e_ES_bg)/(gs[gs_ID].e_sd*sqrt(1/(double)intSize_net+1/(double)bgSize_net));
			P_score_e = z2p(Z_score_e);
			if(P_score_n>1 || P_score_n < 0)
				P_score_n = 1;
			if(P_score_e>1 || P_score_e < 0)
				P_score_e = 1;
			fprintf(fp1,"%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d/%d\t%d/%d\t%d\n",gs_ID,Z_score_n,P_score_n,Z_score_e,P_score_e,gs[gs_ID].size,gs[gs_ID].size_net,intSize,bgSize,intSize_net,bgSize_net,gs[gs_ID].ov);
		}
	}
	fclose(fp1);
	return 1;
}
// output multi results
// LEGOout_int_mul_file((intG,intGNum,intNum,gsNum,geneset_min,geneset_max,output_ES_file,input_NW_file,input_CS_file,input_GS_file);)
LEGOout_int_mul_bg_file(boolean **intG,boolean *bgG,int *intGNum,int *intGNum_net,int bgSize, int bgSize_net, int intNum, int gsNum,int geneset_min, int geneset_max,char *output_ES_file,char *input_NW_file,char *input_CS_file,char *input_GS_file){
	FILE *fp1,*fp2,*fp3,*fp4;
	int i,mm,gs_ID,gene_ID,gs_size,gs_size_net;
	double avg,sd,wei,Z_score_n,Z_score_e,P_score_n,P_score_e;
	boolean *use_gs;
	struct GS_INFO_BG **gs;
	gs  = (struct GS_INFO_BG **)malloc((intNum+1)*sizeof(struct GS_INFO_BG *));
	for(i=1;i<=intNum;i++){
		gs[i]   = (struct GS_INFO_BG *)malloc((gsNum+1)*sizeof(struct GS_INFO_BG));
	}
	use_gs = (boolean *)malloc((gsNum+1)*sizeof(boolean));
	for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
		use_gs[gs_ID]=0;
	}
	if ((fp1 = fopen(output_ES_file,"w")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",output_ES_file, strerror(errno));
		return -1;
	}
	if ((fp2 = fopen(input_NW_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_NW_file, strerror(errno));
		return -1;
	}
	if ((fp3 = fopen(input_CS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_CS_file, strerror(errno));
		return -1;
	}
	if ((fp4 = fopen(input_GS_file,"r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",input_GS_file, strerror(errno));
		return -1;
	}
	// read in GS info file 
	// e.g 1       1       52      0.137187        5.216251	
	while(fscanf(fp4,"%d\t%d\t%d\t%d\t%lf\t%lf\n",&mm,&gs_ID,&gs_size,&gs_size_net,&avg,&sd)!=EOF){
		if(gs_size < geneset_min || gs_size > geneset_max){ // size selection
			continue;
		}
		use_gs[gs_ID] = 1; // use gs
		for(i=1;i<=intNum;i++){
			if(mm==1){ // n-LEGO
				gs[i][gs_ID].size = gs_size;
				gs[i][gs_ID].size_net = gs_size_net;
				gs[i][gs_ID].n_mean = avg;
				gs[i][gs_ID].n_sd = sd;
				gs[i][gs_ID].n_ES = 0;
				gs[i][gs_ID].n_ES_bg = 0;
				gs[i][gs_ID].ov = 0;
			}
			if(mm==2){
				gs[i][gs_ID].e_mean = avg;
				gs[i][gs_ID].e_sd = sd;
				gs[i][gs_ID].e_ES = 0;
				gs[i][gs_ID].e_ES_bg = 0;
			}
		}
	}		
	fclose(fp4);
	printf("finish read in %s\n",input_GS_file);
	// read in NW 
	// 786	1	83.742686
	while(fscanf(fp2,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			for(i=1;i<=intNum;i++){
				if(intG[i][gene_ID]==1){ // belong to int
					gs[i][gs_ID].n_ES += wei;						
					gs[i][gs_ID].ov ++;
				}
				if(bgG[gene_ID]==1){
					gs[i][gs_ID].n_ES_bg += wei;
				}
			}
		}
	}
	fclose(fp2);
	printf("finish read in %s\n",input_NW_file);
	// read in CS
	while(fscanf(fp3,"%d\t%d\t%lf\n",&gene_ID,&gs_ID,&wei)!=EOF){
		if(use_gs[gs_ID]){ // use gs
			for(i=1;i<=intNum;i++){
				if(intG[i][gene_ID]==1){ // belong to int
					gs[i][gs_ID].e_ES += wei;						
				}
				if(bgG[gene_ID]==1){
					gs[i][gs_ID].e_ES_bg += wei;
				}
			}	
		}
	}
	fclose(fp3);
	printf("finish read in %s\n",input_CS_file);
	// output
	double scale_val;
	fprintf(fp1,"Gene_set\tZ_score(n-LEGO)\tP_value(n-LEGO)\tZ_score(e-LEGO)\tP_value(e-LEGO)\tgsSize\tgsSize_net\tintSize\tintSize_net\toverlapSize\tInterestingGeneListID\n");
	for(i=1;i<=intNum;i++){
		for(gs_ID=1;gs_ID<=gsNum;gs_ID++){ // for each gs
			if(use_gs[gs_ID]){ // use gs
				gs[i][gs_ID].n_ES = gs[i][gs_ID].n_ES/intGNum[i];	
				gs[i][gs_ID].n_ES_bg = gs[i][gs_ID].n_ES_bg/bgSize;	
				scale_val = 1/sqrt(1-(double)intGNum[i]/(double)bgSize);
				Z_score_n = scale_val*(gs[i][gs_ID].n_ES-gs[i][gs_ID].n_ES_bg)/(gs[i][gs_ID].n_sd*sqrt(1/(double)intGNum[i]+1/(double)bgSize));
				P_score_n = z2p(Z_score_n);
				gs[i][gs_ID].e_ES = gs[i][gs_ID].e_ES/intGNum_net[i];	
				gs[i][gs_ID].e_ES_bg = gs[i][gs_ID].e_ES_bg/bgSize_net;	
				scale_val = 1/sqrt(1-(double)intGNum_net[i]/(double)bgSize);
				Z_score_e = scale_val*(gs[i][gs_ID].e_ES-gs[i][gs_ID].e_ES_bg)/(gs[i][gs_ID].e_sd*sqrt(1/(double)intGNum_net[i]+1/(double)bgSize_net));
				P_score_e = z2p(Z_score_e);
				if(P_score_n>1 || P_score_n <0)
					P_score_n = 1;
				if(P_score_e>1 || P_score_e < 0)
					P_score_e = 1;
				fprintf(fp1,"%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d/%d\t%d/%d\t%d\t%d\n",gs_ID,Z_score_n,P_score_n,Z_score_e,P_score_e,gs[i][gs_ID].size,gs[i][gs_ID].size_net,intGNum[i],bgSize,intGNum_net[i],bgSize_net,gs[i][gs_ID].ov,i);
			}
		}
	}
	fclose(fp1);
	return 1;
}
#endif
