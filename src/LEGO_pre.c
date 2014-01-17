/* Created by Sherry 
 * Network start with Node:1
 * self-loop will be calculated
 * */

#include "io.h"
#include "graph_basic.h"

int main (int argc, char *argv[]) {	
	char input_net_file[MAX_FILE_NAME_LENGTH];
	char input_gs_file[MAX_FILE_NAME_LENGTH];
	char output_NW_file[MAX_FILE_NAME_LENGTH]; // <output node weight file>
	char output_CS_file[MAX_FILE_NAME_LENGTH]; // <output contribution score file>
	char output_GS_file[MAX_FILE_NAME_LENGTH]; // <output gene set info file>
	double **node_wei; // modified node weight, first column no use
	boolean **node_gs; // matrix indicate whether genes in the gs, start from 1-row, 1-col
	int *gsSize,*gsSize_net,*c1,*c2;
	double *ES,*v1,*v2,bgNE; // edge weight sum
	int i,nodeNum,gsNum;
	struct LGraph lg;
// argc	
	if (argc <= 2) {
		fprintf(stderr, "Usage: %s <network file> <gene set file> <network bg NE (optional, default:0.25)> \n",argv[0]); 
		exit(0); 
	} 
// argc --> par	
	strcpy(input_net_file, argv[1]); // input network file name
	strcpy(input_gs_file, argv[2]); // input gene set file name
	strcpy(output_NW_file,input_gs_file); // output node weight file
	strcat(output_NW_file,"_NW"); // output node weight file
	strcpy(output_GS_file,input_gs_file); // output contribution score file
	strcat(output_GS_file,"_GS"); // output contribution score file
	strcpy(output_CS_file,input_gs_file); // output gene set info file
	strcat(output_CS_file,"_CS"); // output gene set info file
	bgNE = 0.25;
	if(argc >= 3){bgNE=atof(argv[3]);}
	// read in network file, create into LG format
	GetGS_info(input_gs_file,&nodeNum,&gsNum); // nodeNum: maximun number of nodes 
	CreateLG_file(input_net_file,&lg,&nodeNum);	 // create LGraph from file, nodeNum: max(gs_ID,net_ID)
	printf("Real vertex number (max ID index) %d;edge number %d;exist vertex number %d\n",lg.vexnum,lg.arcnum,nodeNum); // test
	// create data	
	node_wei = (double **)malloc((nodeNum+1)*sizeof(double*)); // store node weight for every node in net
	node_gs  = (boolean **)malloc((nodeNum+1)*sizeof(boolean*)); // store gs info for every node in net
	for(i=0;i<=nodeNum;i++){
		node_gs[i]  = (boolean *)malloc((gsNum+1)*sizeof(boolean));
		node_wei[i] = (double *)malloc((gsNum+1)*sizeof(double));
	}
	gsSize = (int *)malloc((gsNum+1)*sizeof(int)); 
	gsSize_net = (int *)malloc((gsNum+1)*sizeof(int)); 
// read in gene set file e.g 1(gene id)   1(gene set id)
	ReadGS_file(input_gs_file,nodeNum,lg.vexnum,&node_gs,&gsSize,&gsSize_net,gsNum);
// get node weight
	GetNodeWei(lg,node_gs,gsSize_net,gsNum,&node_wei,bgNE);
	Out_NW_CS_file(output_NW_file,output_GS_file,output_CS_file,lg,node_wei,gsNum,node_gs,gsSize,gsSize_net,nodeNum);
	return(0);
}
