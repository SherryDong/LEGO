/* Created by Sherry 
 * Network start with Node:1
 * self-loop will be calculated
 * */

#include "io.h"
#include "graph_basic.h"

int main (int argc, char *argv[]) {	
	char input_gs_file[MAX_FILE_NAME_LENGTH]; // <input gs file>
	char input_NW_file[MAX_FILE_NAME_LENGTH]; // <input node weight file>
	char input_CS_file[MAX_FILE_NAME_LENGTH]; // <input contribution score file>
	char input_GS_file[MAX_FILE_NAME_LENGTH]; // <input gene set info file>
	char input_int_file[MAX_FILE_NAME_LENGTH]; // <input interesting gene list file>
	char input_bg_file[MAX_FILE_NAME_LENGTH]; // <input bg gene list file>
	char output_ES_file[MAX_FILE_NAME_LENGTH];
	boolean *intG,*bgG;
	int i,nodeNum,gsNum,intSize,intSize_net,bgSize,bgSize_net;
	int geneset_min, geneset_max;
// argc	
	if (argc < 8) {
		fprintf(stderr, "Usage: %s <gene set file> <node weight file> <contribution score file> <gene set info file> <interesting gene list file> <minimun gene set size> <maximum gene set size> <backgorund file>\n",argv[0]); 
		exit(0); 
	} 
// argc --> par	
	strcpy(input_gs_file, argv[1]); // input gs file name
	strcpy(input_NW_file, argv[2]); // input NW file name
	strcpy(input_CS_file, argv[3]); // input CS file name
	strcpy(input_GS_file, argv[4]); // input GS file name
	strcpy(input_int_file, argv[5]); // input interesting list file name
	geneset_min = atoi(argv[6]);	
	geneset_max = atoi(argv[7]);	
	strcpy(input_bg_file, argv[8]); // input background file name
	strcpy(output_ES_file,input_int_file); // output file
	strcat(output_ES_file,".out"); // output modified network file
// read in gs
	GetGS_info(input_gs_file,&nodeNum,&gsNum); 
// read in interesting gene list and calculate ES
	intG = (boolean *)malloc((MAX_IG_SIZE+1)*sizeof(boolean));
	bgG = (boolean *)malloc((MAX_NODE_SIZE+1)*sizeof(boolean));
	for(i=0;i<=MAX_IG_SIZE;i++){
		intG[i] = 0;
	}
	for(i=0;i<=MAX_NODE_SIZE;i++){
		bgG[i] = 0;
	}
	ReadIG_int_file(input_int_file,input_GS_file,&intG,&intSize,&intSize_net);	
	ReadIG_int_file(input_bg_file,input_GS_file,&bgG,&bgSize,&bgSize_net);	
	LEGOout_int_bg_file(intG,bgG,intSize,intSize_net,bgSize,bgSize_net,gsNum,geneset_min,geneset_max,output_ES_file,input_NW_file,input_CS_file,input_GS_file);
	printf("Finish results in %s\n",output_ES_file);
	return(0);
}

