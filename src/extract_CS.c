/* Created by Sherry 
 * Network start with Node:1
 * self-loop will be calculated
 * */
#include "io.h"
// $txt = `$main_dir/extract_CS $each_gs\_id_gs_CS $gs_id_num $use_gs_num $use_gs $gene_id_num $use_input_num $use_input`; ## neighbor: gsid\tgeneid\tweight
int main (int argc, char *argv[]) {	
	char file_gs_mod[MAX_FILE_NAME_LENGTH];
	int gs_id_num,gs_id_length,gene_id_num,gene_id_length,d1,d2,each,i;
	int start;
	double val;
	FILE  *fp;
	//
	strcpy(file_gs_mod, argv[1]); // id_gs_CS
	gs_id_num = atoi(argv[2]); // 
	gs_id_length = atoi(argv[3]); 
	boolean *use_gs;
	use_gs   = (boolean *)malloc((gs_id_num+1)*sizeof(boolean));
	for(i=1;i<=gs_id_num;i++){
		use_gs[i]=0;
	}
	for(i=1;i<=gs_id_length;i++){
		each = atoi(argv[3+i]);
		//printf("!%d\n",each);
		use_gs[each] = 1;
	}
	start = 4+gs_id_length;
	gene_id_num = atoi(argv[start]); // 
	gene_id_length = atoi(argv[start+1]); // 
	boolean *use_gene;
	use_gene = (boolean *)malloc((gene_id_num+1)*sizeof(boolean));
	for(i=1;i<=gene_id_num;i++){
		use_gene[i]=0;
	}
	for(i=1;i<=gene_id_length;i++){
		each = atoi(argv[start+1+i]);
		//printf("!%d\n",each);
		use_gene[each] = 1;
	}
	i=0;
	fp = fopen(file_gs_mod, "r");
	while(fscanf(fp,"%d\t%d\t%lf",&d1,&d2,&val)!=EOF){
		if(use_gs[d2]==1 && use_gene[d1]==1){
			printf("%d\t%d\t%f\n",d2,d1,val);
		}				
	}
	fclose(fp);
	return 1;	
}


