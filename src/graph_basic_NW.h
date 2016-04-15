/* Created by Sherry @ 20130422
	structure for node, edge in graph
*/

#ifndef _graph_basic
#define _graph_basic

#include "io.h"

// Create graph from file
int CreateLG_file(char *file, struct LGraph *lg,int *nodeNum) { // create graph for no direction and with weight
	char line[MAX_LINE_SIZE];
	char *token; 
	int i,j;
	double wei;
	int ll;
	int *node_index;
	FILE *fp;
	int v1,v2;
	int max=0;
	struct ArcNode *pi; // edge pointer
	node_index = (int*)malloc((MAX_NODE_SIZE)*sizeof(int));
	if ((fp = fopen(file, "r")) == NULL){
		fprintf(stderr, "Coudn't open file %s; %s\n",file, strerror(errno) );
		return -1;
	}
	printf("Successful read in %s , weighted, undirected\n",file);
	ll=0;
	while(fscanf(fp,"%d\t%d\t%lf",&v1,&v2,&wei)!=EOF){
		max=(max>v1)?max:v1;
		max=(max>v2)?max:v2;
		v1--;v2--;
		if(!(*lg).vertexs[v1].head) (*lg).vertexs[v1].head=NULL; // original pointer to NULL
		if(!(*lg).vertexs[v2].head) (*lg).vertexs[v2].head=NULL; // original pointer to NULL
		node_index[v1] = 1;
		node_index[v2] = 1;

		pi=(struct ArcNode *)malloc(sizeof(struct ArcNode));
		pi->adjvex=v2; 
		pi->weight=wei;
		pi->nextarc=(*lg).vertexs[v1].head;
		(*lg).vertexs[v1].head=pi;

		pi=(struct ArcNode *)malloc(sizeof(struct ArcNode));
		pi->adjvex=v1; 
		pi->weight=wei;
		pi->nextarc=(*lg).vertexs[v2].head;
		(*lg).vertexs[v2].head=pi;

		ll++;
	}
	(*lg).vexnum = max;
	(*lg).arcnum  = ll;
	(*nodeNum) = ((*nodeNum) > max)?(*nodeNum):max;
	return(1);
}
// calculate node weight
int GetNodeWei(struct LGraph lg,boolean **node_gs,int *gsSize_net,int gsNum,double ***node_wei, double bgNE){
	struct ArcNode *pi; // edge pointer
	double wei,expect,d1,d2,d3,d4;
	int i,j,id;
	int n=gsNum+1;
	double **node_deg;
	node_deg = (double **)malloc((lg.vexnum+1)*sizeof(double *));
	pi=(struct ArcNode *)malloc(sizeof(struct ArcNode));
	for(i=0;i<lg.vexnum;i++){ // node in network 
		node_deg[i+1] = (double *)malloc((gsNum+1)*sizeof(double));
		for(j=1;j<=gsNum;j++){ // for each gene set
			node_deg[i+1][j]=0;
			(*node_wei)[i+1][j]=0;
		}
	}
	for(i=0;i<lg.vexnum;i++){ // for each gene in network
		pi=lg.vertexs[i].head;
		while(pi!=NULL){
			wei=pi->weight;
			id=pi->adjvex; // interaction ID
			node_deg[i+1][0] += wei;
			for(j=1;j<=gsNum;j++){ // for each gene set
				if(node_gs[id+1][j]==1){ // interaction ID in gs
					node_deg[i+1][j] += wei; 
				}
			}
			pi=pi->nextarc;
		}
		//printf("node degree for node %d: %f\n",i+1,node_deg[i+1][0]);
		// GANPA
		for(j=1;j<=gsNum;j++){ // for each gene set
			//printf("node degree for node %d gene set %d: %f\n",i+1,j,node_deg[i+1][j]);
			if(node_gs[i+1][j]==1){
				d1 = node_deg[i+1][0]*(gsSize_net[j]-1)/(lg.vexnum-1);
			}else{
				d1 = node_deg[i+1][0]*(gsSize_net[j])/(lg.vexnum-1); // do not belong to gs
			}
			d3 = node_deg[i+1][j];
			if(d3 <= d1 && node_gs[i+1][j]==0){ // wi<0, do not belong to gs
				(*node_wei)[i+1][j]=0;
				//(*node_wei)[i+1][j]=1; // no-weight e-LEGO
			}
			if(d3 <= d1 && node_gs[i+1][j]==1){ // wi<0, belong to gs
				(*node_wei)[i+1][j]=1;
			}
			if(d3 > d1 && node_gs[i+1][j]==0){ // wi>0, do not belong to gs
				//(*node_wei)[i+1][j]=pow(d3-d1+1,bgNE*10); 
				//(*node_wei)[i+1][j]=0; // no-neighbor e-LEGO
				(*node_wei)[i+1][j]=1; // no-weight e-LEGO
			}
			if(d3 > d1 && node_gs[i+1][j]==1){ // wi>0, belong to gs
				//(*node_wei)[i+1][j]=pow(d3-d1+1,bgNE*10);  
				(*node_wei)[i+1][j]=1; // no-weight e-LEGO
			}
			//printf("%d\t%d\t%f\n",i+1,j,(*node_wei)[i+1][j]);
		}
	}
	printf("Finish calculating node weight !\n");
	return(1);
}
// calculate node weight
int GetNodeWei_nw(struct LGraph lg,boolean **node_gs,int *gsSize_net,int gsNum,double ***node_wei, double bgNE){
	struct ArcNode *pi; // edge pointer
	double wei,expect,d1,d2,d3,d4;
	int i,j,id;
	int n=gsNum+1;
	double **node_deg;
	node_deg = (double **)malloc((lg.vexnum+1)*sizeof(double *));
	pi=(struct ArcNode *)malloc(sizeof(struct ArcNode));
	for(i=0;i<lg.vexnum;i++){ // node in network 
		node_deg[i+1] = (double *)malloc((gsNum+1)*sizeof(double));
		for(j=1;j<=gsNum;j++){ // for each gene set
			node_deg[i+1][j]=0;
			(*node_wei)[i+1][j]=0;
		}
	}
	for(i=0;i<lg.vexnum;i++){ // for each gene in network
		pi=lg.vertexs[i].head;
		while(pi!=NULL){
			wei=pi->weight;
			id=pi->adjvex; // interaction ID
			node_deg[i+1][0] += wei;
			for(j=1;j<=gsNum;j++){ // for each gene set
				if(node_gs[id+1][j]==1){ // interaction ID in gs
					node_deg[i+1][j] += wei; 
				}
			}
			pi=pi->nextarc;
		}
		//printf("node degree for node %d: %f\n",i+1,node_deg[i+1][0]);
		// GANPA
		for(j=1;j<=gsNum;j++){ // for each gene set
			//printf("node degree for node %d gene set %d: %f\n",i+1,j,node_deg[i+1][j]);
			if(node_gs[i+1][j]==1){
				d1 = node_deg[i+1][0]*(gsSize_net[j]-1)/(lg.vexnum-1);
			}else{
				d1 = node_deg[i+1][0]*(gsSize_net[j])/(lg.vexnum-1); // do not belong to gs
			}
			d3 = node_deg[i+1][j];
			(*node_wei)[i+1][j]=1; // no-weight e-LEGO
			//printf("%d\t%d\t%f\n",i+1,j,(*node_wei)[i+1][j]);
		}
	}
	printf("Finish calculating node weight !\n");
	return(1);
}
// calculate node weight
int GetNodeWei_nn(struct LGraph lg,boolean **node_gs,int *gsSize_net,int gsNum,double ***node_wei, double bgNE){
	struct ArcNode *pi; // edge pointer
	double wei,expect,d1,d2,d3,d4;
	int i,j,id;
	int n=gsNum+1;
	double **node_deg;
	node_deg = (double **)malloc((lg.vexnum+1)*sizeof(double *));
	pi=(struct ArcNode *)malloc(sizeof(struct ArcNode));
	for(i=0;i<lg.vexnum;i++){ // node in network 
		node_deg[i+1] = (double *)malloc((gsNum+1)*sizeof(double));
		for(j=1;j<=gsNum;j++){ // for each gene set
			node_deg[i+1][j]=0;
			(*node_wei)[i+1][j]=0;
		}
	}
	for(i=0;i<lg.vexnum;i++){ // for each gene in network
		pi=lg.vertexs[i].head;
		while(pi!=NULL){
			wei=pi->weight;
			id=pi->adjvex; // interaction ID
			node_deg[i+1][0] += wei;
			for(j=1;j<=gsNum;j++){ // for each gene set
				if(node_gs[id+1][j]==1){ // interaction ID in gs
					node_deg[i+1][j] += wei; 
				}
			}
			pi=pi->nextarc;
		}
		//printf("node degree for node %d: %f\n",i+1,node_deg[i+1][0]);
		// GANPA
		for(j=1;j<=gsNum;j++){ // for each gene set
			//printf("node degree for node %d gene set %d: %f\n",i+1,j,node_deg[i+1][j]);
			if(node_gs[i+1][j]==1){
				d1 = node_deg[i+1][0]*(gsSize_net[j]-1)/(lg.vexnum-1);
			}else{
				d1 = node_deg[i+1][0]*(gsSize_net[j])/(lg.vexnum-1); // do not belong to gs
			}
			d3 = node_deg[i+1][j];
			if(d3 <= d1 && node_gs[i+1][j]==0){ // wi<0, do not belong to gs
				(*node_wei)[i+1][j]=0;
			}
			if(d3 <= d1 && node_gs[i+1][j]==1){ // wi<0, belong to gs
				(*node_wei)[i+1][j]=1;
			}
			if(d3 > d1 && node_gs[i+1][j]==0){ // wi>0, do not belong to gs
				(*node_wei)[i+1][j]=0; // no-neighbor e-LEGO
			}
			if(d3 > d1 && node_gs[i+1][j]==1){ // wi>0, belong to gs
				(*node_wei)[i+1][j]=pow(d3-d1+1,bgNE*10);  
			}
			//printf("%d\t%d\t%f\n",i+1,j,(*node_wei)[i+1][j]);
		}
	}
	printf("Finish calculating node weight !\n");
	return(1);
}
#endif
