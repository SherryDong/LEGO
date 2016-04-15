/***

iNP.C: An iterative Network Partition (iNP) algorithm
authors: Weidong Tian & Siqi Sun
School of Life Sciences, Institute of Biostatistics
Fudan University
email: weidong.tian@fudan.edu.cn
last modified: 08/05/2011

***/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>

#define MAX_DIR 200
#define BIN_DIR "/home/sherry/project/network/NETGO/LEGO_code_v6.0_RNA/bin//"
#define MATRIX_DIR ""
#define LOWER_BOUND 5
#define MAX_ATTR_NAME 400

typedef struct _Vertex *VertexPtr;    // pointer to vertex structure 
typedef struct _Vertex {
	char name[MAX_ATTR_NAME];
	int index; 
} VertexStr;
VertexPtr vertices;

typedef struct _Edge *EdgePtr;    // pointer to edge structure 
typedef struct _Edge {
	int first;
	int second;
} EdgeStr;

typedef struct _Module *ModulePtr; // pointer to network module structure
typedef struct _Module {
	int num_nodes; 
	int *nodes; 
} ModuleStr;  

char **edges;       // global variables of network edges
ModulePtr final_modules;   // pointer of final partitioned modules 
int num_final_modules, *pnum_final_modules; // number of final partitioned modules 
int tt; 

/* function to partition networks */
double partition_network(char input_matrix_file[], int *pnum_modules, ModulePtr modules, char prog[], char prog_option); 

/* function for recursive partitioning */
void iRun (ModulePtr module, int num_vertices, char prog[], char input_matrix_name[], char prog_option); 

/* function to add new modules in the recursive partition */ 
void add_modules (ModulePtr module);

/* function to process the input network file */ 
int process_input (char input_file[], char input_index_file[], char input_matrix_file[], int *pnum_vertices);

/* function to print final partitioned modules */ 
int print_modules(char final_module_file[], double M); 
int get_index (char name[MAX_ATTR_NAME], int num_vertices);

/* function to read input network matrix */ 
int read_input_matrix (char input_matrix_file[], int num_vertices, int *degrees);
int cal_M_th (double M, int node_size, int num_vertices, char mode, double tmp_M, double th1_up, double th1_low, double th2_up, double th2_low); 

/* function to get random numbers */ 
int GetRand (int M); 

/* function to compute modularity */ 
double calculate_modularity (int num_modules, ModulePtr modules, int *degrees, int total_edges);
char choose_mode (double M); 

/* function to read modules */ 
double read_modules(char answer_file[], int *pnum_modules, ModulePtr modules);


int process_file (char input_file[]);

/* function to compute Jaccard accuracy */ 
double jaccard (char answer_file[], char input_file[]); 

/* function to print random graphs */ 
int print_rand_network(char tmp_matrix_file[], EdgePtr links, int num_links); 

/* function to partition random graphs */ 
double partition_rand_network(char input_matrix_file[]); 
int random_graph_M (char tmp_matrix_file[], ModulePtr module, char prog[], char prog_option, double M); 

/* function to shuffle edges in order to generate random graphs */ 
void Shuffle (EdgePtr links, char **tmp_edges, int num_links);

int main (int argc, char *argv[]) {
	char input_file[MAX_DIR] = MATRIX_DIR; 
	char input_index_file[MAX_DIR] = MATRIX_DIR; 
	char input_matrix_file[MAX_DIR] = MATRIX_DIR; 
	char input_matrix_name[MAX_DIR];
	char tmp_module_file[MAX_DIR] = MATRIX_DIR; 
	char final_module_file[MAX_DIR] = MATRIX_DIR; 
	char output_file[MAX_DIR] = MATRIX_DIR;
	char log_file[MAX_DIR] = MATRIX_DIR;
	char answer_file[MAX_DIR] = MATRIX_DIR;
	char prog[MAX_DIR]=BIN_DIR; 
	char mode, prog_option; 
	int i,j, k, l, num_edges, *pnum_edges, num_modules, *pnum_modules, num_answer_modules, *pnum_answer_modules; 
	ModulePtr modules, module, answer_modules, answer_module; 
	char tmp_file[MAX_DIR] = MATRIX_DIR; 
	double M, final_M; 
	int num_vertices, *pnum_vertices, num_lines; 
	int *degrees, total_edges; 
	int rand_num; 
	FILE *fp, *log; 
	time_t t_0, t_1, t, t_last;
	clock_t clo_0, clo_1, clo_last; 
	int elapTicks; 
	double elapMilli, elapSeconds, elapMinutes;  
	char answer_option; 
	double acc_0, acc_1; 
	double th1_up, th1_low, th2_up, th2_low; 
	char suffix[MAX_ATTR_NAME]; 
	char rand_graph_option; 

	if (argc < 3) { 
		fprintf(stderr, "Too few arguments!\n");
		fprintf(stderr, "Usage:\n");
		fprintf(stderr, "<Program Name> input_matrix (tab delimitted (string or index)) program_option (1: MSG.out; 2: SCNewman; 3: testQcut) answer_option (0: no answer; 1: answer) answer_file\n"); 
		exit(EXIT_FAILURE);
	}

	strcat(input_file, argv[1]); 
	strcpy (input_matrix_name, argv[1]); 
	prog_option = atoi(argv[2]); 
	answer_option = atoi(argv[3]); 
	if (answer_option) { 
		if (argc < 4) { 
			fprintf(stderr, "Missing answer file\n"); 
			fprintf(stderr, "Usage:\n");
			fprintf(stderr, "<Program Name> input_matrix (tab delimitted (string or index)) program_option (1: MSG.out; 2: SCNewman; 3: testQcut) answer_option (0: no answer; 1: answer) answer_file\n");
			exit(EXIT_FAILURE);
		}
		strcat(answer_file, argv[4]); 
	}

	strcat(input_index_file, argv[1]); 
	strcat(input_index_file, ".index"); 
	strcat(input_matrix_file, argv[1]); 
	strcat(input_matrix_file, ".matrix"); 

	strcat(output_file, argv[1]); 
	strcat(output_file, "_result_"); 
	strcat(log_file, argv[1]); 
	strcat(log_file, "_"); 

	strcat(final_module_file, argv[1]); 
	strcat(final_module_file, "_result_"); 

	if (prog_option == 1) { 
		strcat (prog, "MSG.out"); 
		strcat (input_matrix_name, "_MSG.out"); 
		strcat(final_module_file, "iNP_MSG"); 
		strcat(output_file, "MSG"); 
		strcat(log_file, "iNP_MSG.log"); 
	} else if (prog_option == 2) { 
		strcat (prog, "do_SCNewman"); 
		strcat (input_matrix_name, "_SCNewman");
		strcat(final_module_file, "iNP_SCNewman"); 
		strcat(output_file, "SCNewman"); 
		strcat(log_file, "iNP_SCNewman.log"); 
	} else if (prog_option == 3) { 
		strcat (prog, "qcut.pl");
		strcat (input_matrix_name, "_qcut"); 
		strcat(final_module_file, "iNP_Qcut"); 
		strcat(output_file, "Qcut"); 
		strcat(log_file, "iNP_Qcut.log"); 
	} else { 
		fprintf(stderr, "<Program Name> input_matrix (tab delimitted (string or index)) mode (0, 1, 2) program_option (1: MSG.out)\n"); 
		fprintf(stderr, "program_option currently allows only 1 - 3 (1: MSG, 2: SC, 3: Qcut)\n"); 
		exit(EXIT_FAILURE); 
	}
	if ((log = fopen(log_file, "w")) == NULL) {
		fprintf(stderr, "\nFailure to write file %s in read mode\n", log_file);
		fflush(NULL);
		return(-1);
	}

	num_vertices = 0; 
	pnum_vertices = &num_vertices;
	num_lines =	process_input(input_file, input_index_file, input_matrix_file, pnum_vertices); 
	degrees = (int *) calloc(num_vertices, sizeof(int))-1;
	total_edges = read_input_matrix(input_matrix_file, num_vertices, degrees); 

	num_modules = 0; 
	pnum_modules = &num_modules; 
	modules = (ModulePtr) malloc(num_vertices * sizeof(ModuleStr)) - 1; 

	time (&t_0); 
	M = partition_network(input_matrix_file, pnum_modules, modules, prog, prog_option); 
	time (&t_1); 

	if ((fp = fopen(output_file, "w")) == NULL) {
		fprintf(stderr, "\nFailure to write file %s in read mode\n", output_file);
		fflush(NULL);
		return(-1);
	}

/* print networks partitioned without iterations */ 
	for (i = 1; i <= num_modules; i ++) { 
		module = modules + i; 
		for (j = 1; j <= module->num_nodes; j ++) { 
			fprintf(fp, "%s\t", (vertices + module->nodes[j])->name); 
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "%f\n", M); 
	fclose (fp); 

	final_modules = (ModulePtr) malloc(num_vertices * sizeof(ModuleStr)) - 1; 
	num_final_modules = 0; 
	pnum_final_modules = &num_final_modules; 

	//printf ("%s %f\n", output_file, M); 
	rand_num = (int) rand(); 
	for (i = 1; i <= num_modules; i ++) { 
		module = modules + i; 
		iRun(module, num_vertices, prog, input_matrix_name, prog_option); 
	}
	if (num_final_modules) { 
		final_M = calculate_modularity(num_final_modules, final_modules, degrees, total_edges); 
		/* print networks partitioned with iterations */ 
		print_modules(final_module_file, final_M); 
	}
	time (&t_last); 

	fprintf (log, "NP\tnum_modules\t%d\tmodularity\t%f\ttime\t%f\t", num_modules, M, difftime(t_1, t_0)); 
	if (answer_option) { 
		acc_0 = jaccard(answer_file, output_file); 
		fprintf (log, "acc\t%f\n", acc_0); 
	} else { 
		fprintf (log, "\n"); 
	}
	fprintf (log, "iNP\tnum_modules\t%d\tmodularity\t%f\ttime\t%f\t", num_final_modules, final_M, difftime(t_last, t_0)); 
	if (answer_option) { 
		acc_1 = jaccard(answer_file, final_module_file); 
		fprintf (log, "acc\t%f\n", acc_1); 
	} else { 
		fprintf (log, "\n"); 
	}
	fclose (log); 
	
	remove(input_matrix_file); 
	remove(input_index_file); 

	return 0;
}

/* function to compute the jaccard accuracy of comparing a partitioned network with a defined network */ 
double jaccard (char answer_file[], char input_file[]) { 
	double sim, max, sum1, sum2, acc;
	int overlap, num_lines, num_nodes;
	char *tmp_nodes;
	int i,j, k, l, num_edges, *pnum_edges, num_modules, *pnum_modules, num_answer_modules, *pnum_answer_modules;
	ModulePtr modules, module, answer_modules, answer_module;

	num_answer_modules = 0;
	pnum_answer_modules = &num_answer_modules;
	num_lines = process_file(answer_file);
	answer_modules = (ModulePtr) malloc(num_lines * sizeof(ModuleStr)) - 1;

	num_modules = 0;
	pnum_modules = &num_modules;
	num_lines = process_file(input_file);
	modules = (ModulePtr) malloc(num_lines * sizeof(ModuleStr)) - 1;

	read_modules(answer_file, pnum_answer_modules, answer_modules);
	read_modules(input_file, pnum_modules, modules);

	num_nodes = 0;
	for (i = 1; i <= num_modules; i ++) {
		module = modules + i;
		for (j = 1; j <= module->num_nodes; j ++) {
			if (num_nodes < module->nodes[j]) {
				num_nodes = module->nodes[j];
			}
		}
	}
	sum1 = 0; sum2 = 0;
	for (i = 1; i <= num_modules; i ++) {
		module = modules + i;
		tmp_nodes = (char *) calloc(num_nodes, sizeof(char)) - 1;
		for (j = 1; j <= module->num_nodes; j ++) {
			tmp_nodes[module->nodes[j]] = 1;
		}
		max = 0;
		if (module->num_nodes == 0) { continue;}
		for (l = 1; l <= num_answer_modules; l ++) {
			answer_module = answer_modules + l;
			overlap = 0;
			for (j = 1; j <= answer_module->num_nodes; j ++) {
				if (tmp_nodes[answer_module->nodes[j]]) {
					overlap ++;
				}
			}
			sim = (double) overlap / (module->num_nodes + answer_module->num_nodes - overlap);
			if (max < sim) {
				max = sim;
			}
		}
		sum1 += max * module->num_nodes;
		sum2 += module->num_nodes;
		free(tmp_nodes + 1);
	}
	acc = (double) sum1 / sum2;
	for (i = 1; i <= num_modules; i ++) {
		module = modules + i;
		if (module->num_nodes == 0) { continue;}
		free (module->nodes + 1); 
	}
	for (i = 1; i <= num_answer_modules; i ++) {
		module = answer_modules + i;
		if (module->num_nodes == 0) { continue;}
		free (module->nodes + 1); 
	}
	free(answer_modules + 1); 
	free(modules + 1); 
	return acc;
}

int process_file (char input_file[]) {
	FILE *fp, *out_index, *out_matrix;
	int line_number;
	char *token1, *token2, *token3;
	char line[138327];  // length 2000 
	double a;
	char first_line;
	char mark[40];
	int i, j, max_id, num_lines, index, num_names;
	char name1[MAX_ATTR_NAME];
	char name2[MAX_ATTR_NAME];

	if ((fp = fopen(input_file, "r")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", input_file);
		fflush(NULL);
		return(-1);
	}
	i = 0; j = 0; num_lines = 0; index = 0;
	while (fgets(line, sizeof(line), fp) != NULL) {
		num_lines ++;
	}
	fclose(fp);
	return num_lines;
}


double read_modules(char result_file[], int *pnum_modules, ModulePtr modules) {
	FILE *fp, *out_index, *out_matrix;
	int line_number;
	char *token1, *token2, *token3;
	char line[138327];  
	double M;
	char first_line;
	int i, j, max_id, num_lines, index, num_names;
	char command[MAX_ATTR_NAME];
	char mark[40];
	ModulePtr module;

	if ((fp = fopen(result_file, "r")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", result_file);
		fflush(NULL);
		return(-1);
	}

	index = 0;
	strcpy(mark, "0");
	M = 0;
	while (fgets(line, sizeof(line), fp) != NULL) {
		token1 = strtok(line, "\t");
		if (token1 == NULL) {continue;}
		if ((strncmp(token1, mark, 1)) == 0) {
			M = atof(token1);
			break;
		}
		index ++;
		(*pnum_modules) = index;
		module = modules + index;
		module->num_nodes = 1;
		module->nodes = (int *) malloc(sizeof(int)) - 1;
		module->nodes[1] = atoi(token1);
		while ((token1 = strtok(NULL, "\t\n")) != NULL ) {
			module->num_nodes ++;
			module->nodes = (int *) realloc(module->nodes + 1, module->num_nodes * sizeof(int)) - 1;
			module->nodes[module->num_nodes] = atoi(token1);
		}
	}
	fclose (fp);
	if (index) {
		modules = (ModulePtr) realloc(modules + 1, index * sizeof(ModuleStr)) - 1;
	}
	return M;
}

int read_input_matrix (char input_matrix_file[], int num_vertices, int *degrees) {
	FILE *fp, *out_index, *out_matrix;
	int line_number;
	char *token1, *token2, *token3;
	char line[138327];  
	double a;  
	char first_line;
	char mark[40]; 
	int i, j, max_id, num_lines, index, num_names; 
	char name1[MAX_ATTR_NAME];
	char name2[MAX_ATTR_NAME];

	if ((fp = fopen(input_matrix_file, "r")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", input_matrix_file);
		fflush(NULL);
		return(-1);
	}

	edges = (char **) calloc(num_vertices, sizeof(char *)) - 1; 
	for (i = 1; i <= num_vertices; i ++) { 
		edges[i]= (char *) calloc(num_vertices, sizeof(char)) - 1; 
	}
	i = 0; j = 0; num_lines = 0; index = 0; 
	while (fgets(line, sizeof(line), fp) != NULL) { 
		token1 = strtok(line, "\t"); 
		if (token1 == NULL) {continue;}
		token2 = strtok(NULL, "\t"); 
		if (token2 == NULL) {continue;}
		i = atoi(token1); j = atoi(token2); 
		if (i == j) {continue;}
		num_lines ++; 
		edges[i][j] = 1;
		edges[j][i] = 1;
	}
	fclose (fp);
	
	num_lines = 0; 
	for (i = 1; i <= num_vertices; i ++) { 
		for (j = 1; j <= num_vertices; j ++) { 
			if (edges[i][j]) { 
				degrees[i] ++;
				num_lines ++; 
			}
		}
	}

	return num_lines/2; 
}

int process_input (char input_file[], char input_index_file[], char input_matrix_file[], int *pnum_vertices) {
	FILE *fp, *out_index, *out_matrix;
	int line_number;
	char *token1, *token2, *token3;
	char line[138327];  
	double a;  
	char first_line;
	char mark[40]; 
	int i, j, max_id, num_lines, index, num_names; 
	char name1[MAX_ATTR_NAME];
	char name2[MAX_ATTR_NAME];

	if ((fp = fopen(input_file, "r")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", input_file);
		fflush(NULL);
		return(-1);
	}
	if ((out_index = fopen(input_index_file, "w")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", input_index_file);
		fflush(NULL);
		return(-1);
	}
	if ((out_matrix = fopen(input_matrix_file, "w")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", input_matrix_file);
		fflush(NULL);
		return(-1);
	}
	i = 0; j = 0; num_lines = 0; index = 0; 
	while (fgets(line, sizeof(line), fp) != NULL) { 
		token1 = strtok(line, "\t"); 
		if (token1 == NULL) {continue;}
		strncpy(name1, token1, MAX_ATTR_NAME);
		token2 = strtok(NULL, "\t"); 
		if (token2 == NULL) {continue;}
		strncpy(name2, token2, MAX_ATTR_NAME);
		if (strcmp(name1, name2) == 0) { continue;}
		token3 = strtok(NULL, "\t\n"); 
		if (token3 == NULL) {continue;}
		num_lines ++; 
		if (num_lines == 1) { 
			vertices = (VertexPtr) malloc(2 * sizeof(VertexStr)) - 1; 
			index ++; 
			i = index; 
			(vertices + index)->index = index; 
			strcpy((vertices + index)->name, name1); 
			fprintf(out_index, "%d\t%s\n", index, name1); 
			index ++; 
			j = index; 
			(vertices + index)->index = index; 
			strcpy((vertices + index)->name, name2); 
			fprintf(out_index, "%d\t%s\n", index, name2); 
			fprintf(out_matrix, "%d\t%d\t%s\n", i, j, token3); 
			(*pnum_vertices) = index; 
		} else { 
			i = get_index(name1, index); 
			if (i > index) { index = i; fprintf(out_index, "%d\t%s\n", index, name1);}
			j = get_index(name2, index); 
			if (j > index) { index = j; fprintf(out_index, "%d\t%s\n", index, name2);}
			fprintf(out_matrix, "%d\t%d\t%s\n", i, j, token3); 
			(*pnum_vertices) = index; 
		}
	}
	fclose(fp); 
	fclose(out_index); 
	fclose(out_matrix); 
	return num_lines; 
}

int get_index (char name[MAX_ATTR_NAME], int num_vertices) { 
	int i, j, index; 
	index = 0; 
	for (i = 1; i <= num_vertices; i ++) { 
		if (strcmp((vertices + i)->name, name) == 0) { 
			index = i; 
			break;
		}
	}
	if (index == 0) { 
		num_vertices ++; 
		vertices = (VertexPtr) realloc(vertices + 1, num_vertices * sizeof(VertexStr)) - 1; 
		(vertices + num_vertices)->index = num_vertices; 
		strcpy((vertices + num_vertices)->name, name);
		index = num_vertices; 
	}
	return index; 
}

double partition_network(char input_matrix_file[], int *pnum_modules, ModulePtr modules, char prog[], char prog_option) { 
	FILE *fp, *out_index, *out_matrix;
	int line_number;
	char *token1, *token2, *token3;
	char line[138327];  // length 2000 
	double M; 
	char first_line;
	int i, j, max_id, num_lines, index, num_names; 
	char command[MAX_ATTR_NAME]; 
	char result_file[MAX_DIR] = MATRIX_DIR; 
	char mark[40]; 
	ModulePtr module; 

	strcpy (command, prog); 
	strcat (command, " "); 
	strcat (command, input_matrix_file); 
	strcat (command,  " 1>/dev/null"); 

	if (prog_option == 1) { 
		strcat (result_file, "MSG_result_of_"); 
	} else if (prog_option == 2) { 
		strcat (result_file, "spectral_result_of_"); 
	} else if (prog_option == 3) { 
		strcat (result_file, "qcut_result_of_"); 
	}
	strcat (result_file, input_matrix_file); 
	//printf ("%s\n", command); 
	system (command); 
	
	if ((fp = fopen(result_file, "r")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", result_file);
		fflush(NULL);
		return(-1);
	}

	index = 0; 
	strcpy(mark, "0");
	while (fgets(line, sizeof(line), fp) != NULL) { 
		token1 = strtok(line, "\t"); 
		if (token1 == NULL) {continue;}
		if ((strncmp(token1, mark, 1)) == 0) { 
			M = atof(token1); 
			break;
		}
		if (atoi(token1) < 1) { continue;}
		index ++; 
		(*pnum_modules) = index; 
		module = modules + index; 
		module->num_nodes = 1; 
		module->nodes = (int *) malloc(sizeof(int)) - 1; 
		module->nodes[1] = atoi(token1); 
		while ((token1 = strtok(NULL, "\t\n")) != NULL ) {
			module->num_nodes ++; 
			module->nodes = (int *) realloc(module->nodes + 1, module->num_nodes * sizeof(int)) - 1; 
			module->nodes[module->num_nodes] = atoi(token1); 
		}
	}
	fclose (fp); 
	remove(result_file); 
	if (index) { 
		modules = (ModulePtr) realloc(modules + 1, index * sizeof(ModuleStr)) - 1; 
	}
	return M;
}

int print_modules(char final_module_file[], double M) { 
	FILE *fp; 
	int i, j; 
	int vertex_1, vertex_2; 
	ModulePtr module; 
	
	if ((fp = fopen(final_module_file, "w")) == NULL) {
		fprintf(stderr, "\nFailure to write file %s in read mode\n", final_module_file);
		fflush(NULL);
		return(-1);
	}

	for (i = 1; i <= num_final_modules; i ++) { 
		module = final_modules + i; 
		if (module->num_nodes) { 
			for (j = 1; j <= module->num_nodes; j ++) { 
				fprintf(fp, "%s\t", (vertices + module->nodes[j])->name); 
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "%f\n", M); 
	fclose (fp); 
	return 0; 
}

int print_matrix(ModulePtr module, char tmp_matrix_file[]) { 
	FILE *fp; 
	int i, j; 
	int vertex_1, vertex_2; 
	int is_isolated, is_singular; 

	if ((fp = fopen(tmp_matrix_file, "w")) == NULL) {
		fprintf(stderr, "\nFailure to write file %s in read mode\n", tmp_matrix_file);
		fflush(NULL);
		return(-1);
	}

	for (i = 1; i < module->num_nodes; i ++) { 
		vertex_1 = module->nodes[i]; 
		for (j = i + 1; j <= module->num_nodes; j ++) { 
			vertex_2 = module->nodes[j]; 
			if (edges[vertex_1][vertex_2]) { 
				fprintf(fp, "%d\t%d\t1\n", vertex_1, vertex_2); 
			}
		}
	}
	fclose (fp); 
	return 0; 
}

void add_modules (ModulePtr module) { 
	ModulePtr tmp_module; 
	int i; 

	(*pnum_final_modules) ++; 
	tmp_module = final_modules + (*pnum_final_modules); 
	tmp_module->num_nodes = module->num_nodes; 
	tmp_module->nodes = (int *) malloc (module->num_nodes * sizeof(int)) - 1; 
	for (i = 1; i <= module->num_nodes; i ++) { 
		tmp_module->nodes[i] = module->nodes[i]; 
	}
	return; 
}

void iRun (ModulePtr module, int num_vertices, char prog[], char input_matrix_name[], char prog_option) { 
	double tmp_M; 
	char tmp_matrix_file[MAX_DIR]= MATRIX_DIR; 
	char tmp_module_file[MAX_DIR]= MATRIX_DIR; 
	char suffix[MAX_ATTR_NAME]; 
	ModulePtr tmp_modules, tmp_module;
	int num_tmp_modules, *pnum_tmp_modules, *pnum_final_modules; 
	int i; 
	double M_th; 
	int pid;
	char test_result, is_singular; 

	if (module->num_nodes <= LOWER_BOUND) { 
		add_modules(module); 
		return; 
	}
	strcat (tmp_matrix_file, input_matrix_name); 
	is_singular = print_matrix(module, tmp_matrix_file); 
	if (is_singular) { 
		add_modules(module); 
		return;
	}
	num_tmp_modules = 0; 
	pnum_tmp_modules = &num_tmp_modules; 
	tmp_modules = (ModulePtr) malloc(module->num_nodes * sizeof(ModuleStr)) - 1; 
	tmp_M = partition_network(tmp_matrix_file, pnum_tmp_modules, tmp_modules, prog, prog_option); 
	remove(tmp_matrix_file); 
	test_result = 0; 
	if (tmp_M >= 0.3) { 
		test_result = random_graph_M (tmp_matrix_file, module, prog, prog_option, tmp_M);
	}
	if (!test_result) { 
		for (i = 1; i <= num_tmp_modules; i ++) { 
			tmp_module = tmp_modules + i; 
			free (tmp_module->nodes + 1); 
		}
		free (tmp_modules + 1); 
		add_modules(module); 
		return;
	}
	for (i = 1; i <= num_tmp_modules; i ++) { 
		tmp_module = tmp_modules + i; 
		iRun(tmp_module, num_vertices, prog, input_matrix_name, prog_option);
	}
	return; 
}

double calculate_modularity (int num_modules, ModulePtr modules, int *degrees, int total_edges) {
	int i, j, k; 
	int *sum_degree, *sum_edges; 
	ModulePtr module; 
	int vertex_1, vertex_2; 
	double M, r1, r2; 

	sum_degree = (int *)calloc(num_modules, sizeof(int)) - 1; 
	sum_edges = (int *)calloc(num_modules, sizeof(int)) - 1; 
	M = 0.0; 
	for (i = 1; i <= num_modules; i ++) { 
		module = modules + i; 
		for (j = 1; j <= module->num_nodes; j ++) { 
			vertex_1 = module->nodes[j]; 
			sum_degree[i] += degrees[vertex_1]; 
			if ( j == module->num_nodes) { break;}
			for (k = j + 1; k <= module->num_nodes; k ++) { 
				vertex_2 = module->nodes[k]; 
				if (edges[vertex_1][vertex_2]) { 
					sum_edges[i] ++; 
				}
			}
		}
		//printf ("module%d sum_degree %d sum_edges %d total_edges %d\n", i, sum_degree[i], sum_edges[i], total_edges); 
		r1 = (double) sum_edges[i] / total_edges; 
		r2 = (double) sum_degree[i] / (2 * total_edges); 
		r2 = r2 * r2; 
		M += r1 - r2; 
	}
	free(sum_degree + 1); 
	free(sum_edges + 1); 
	return M; 
}

int GetRand (int M) { 
	int rand_num; 
	double d; 

	d = (double) rand() / ((double) RAND_MAX + 1); 
	rand_num = (int) (d * M); 

	return (rand_num + 1); 
}

int random_graph_M (char tmp_matrix_file[], ModulePtr module, char prog[], char prog_option, double M) { 
	int i, j, k, l, m, n, total_edges; 
	EdgePtr links; 
	double tmp_M, z, sd, avr, avr2; 
	int num_tmp_modules, *pnum_tmp_modules, *pnum_final_modules; 
	int *nodes;
	char **tmp_edges;
	int iter; 
	double r; 

	nodes = (int *) calloc(module->num_nodes, sizeof(int)) - 1; 
	tmp_edges = (char **) calloc(module->num_nodes, sizeof(char*)) - 1; 
	for (i = 1; i <= module->num_nodes; i ++) { 
		nodes[i] = module->nodes[i]; 
		tmp_edges[i] = (char *) calloc(module->num_nodes, sizeof(char)) - 1; 
	}
	n = module->num_nodes * (module->num_nodes - 1) / 2; 
	links = (EdgePtr) calloc(n, sizeof(EdgeStr)) - 1; 

	avr = 0; avr2 = 0; 
	iter = 50; r = 0;  
	for (l = 1; l <= iter; l ++) {
		m = 0; 
		for (i = 1; i <= module->num_nodes - 1; i ++) { 
			for (j = i + 1; j <= module->num_nodes; j ++) { 
				tmp_edges[i][j] = 0; 
				tmp_edges[j][i] = 0; 
				if (edges[module->nodes[i]][module->nodes[j]]) { 
					tmp_edges[i][j] = 1; 
					tmp_edges[j][i] = 1; 
					m ++;
					(links+m)->first = i;
					(links+m)->second = j;
				}
			}
		}
		total_edges = m;
		Shuffle(links, tmp_edges, total_edges);
		print_rand_network(tmp_matrix_file, links, total_edges); 
		tmp_M = partition_rand_network(tmp_matrix_file); 
		r += M - tmp_M; 
		
		if ( l < 10) {
			if (r >= 1.0) {return 1; }
			if (r < 0) {return 0; }
		} else if (l == 10) { 
			if (r >= 1.0) { return 1; }
			if (r <= 0.2) { return 0;}
		}
		avr += tmp_M;
		avr2 += tmp_M * tmp_M; 
		//printf ("%d -> M %f r %f\n", l, tmp_M, r); 
	}
	
	avr /= (double) iter; 
	avr2 /= (double) iter; 
	sd = sqrt(avr2 - avr * avr); 
	if (sd == 0) { 
		if (M > avr) { 
			return 1; 
		} else { 
			return 0;
		}
	} else { 
		z = (M-avr)/sd; 
	}
	//printf ("M %f avr %f z %f\n", M, avr, z); 
	
	free(links+1); 
	free(nodes+1); 
	for (i = 1; i <= module->num_nodes; i ++) { 
		free(tmp_edges[i]+1); 
	}
	free(tmp_edges + 1); 
	if (z < 2) { 
		return 0; 
	} else { 
		return 1; 
	}
}

int print_rand_network(char tmp_matrix_file[], EdgePtr links, int num_links) { 
	FILE *fp; 
	int i, j; 
	int vertex_1, vertex_2; 

	if ((fp = fopen(tmp_matrix_file, "w")) == NULL) {
		fprintf(stderr, "\nFailure to write file %s in read mode\n", tmp_matrix_file);
		fflush(NULL);
		return(-1);
	}

	for (i = 1; i <= num_links; i ++) { 
		fprintf(fp, "%d\t%d\t1\n", (links+i)->first, (links+i)->second); 
	}
	fclose (fp); 
	return 0; 
}

double partition_rand_network(char input_matrix_file[]) { 
	FILE *fp, *out_index, *out_matrix;
	int line_number;
	char *token1, *token2, *token3;
	char line[138327];  // length 2000 
	double M; 
	char first_line;
	int i, j, max_id, num_lines, index, num_names; 
	char command[MAX_ATTR_NAME]; 
	char result_file[MAX_DIR] = MATRIX_DIR; 
	char mark[40]; 
	ModulePtr module; 
	char prog[MAX_DIR]=BIN_DIR; 

	strcat (prog, "MSG.out"); 
	strcpy (command, prog); 
	strcat (command, " "); 
	strcat (command, input_matrix_file); 
	strcat (command,  " 1>/dev/null"); 

	strcat (result_file, "MSG_result_of_"); 
	strcat (result_file, input_matrix_file); 
	//printf ("%s\n", command); 
	system (command); 
	
	if ((fp = fopen(result_file, "r")) == NULL) {
		fprintf(stderr, "\nFailure to open file %s in read mode\n", result_file);
		fflush(NULL);
		return(-1);
	}

	index = 0; 
	strcpy(mark, "0");
	while (fgets(line, sizeof(line), fp) != NULL) { 
		token1 = strtok(line, "\t"); 
		if (token1 == NULL) {continue;}
		if ((strncmp(token1, mark, 1)) == 0) { 
			M = atof(token1); 
			break;
		}
	}
	fclose (fp); 
	remove(result_file); 
	return M;
}

void Shuffle (EdgePtr links, char **tmp_edges, int num_links) {
	int loop, i, j, t1, t2, n1, n2, m1, m2;
	EdgePtr edge1, edge2, temp1, temp2;

	loop = 10 * num_links;
	for(i=0; i<loop; i++) {
		t1 = GetRand(num_links);
		t2 = GetRand(num_links);
		while (t2 == t1) {
			t2 = GetRand(num_links);
		}
		edge1 = links + t1;
		edge2 = links + t2;
		n1 = edge1->first;
		n2 = edge1->second;
		m1 = edge2->first;
		m2 = edge2->second;
		if (tmp_edges[n1][m2] || tmp_edges[m1][n2] || n1 == m2 || n2 == m1) {
			continue;
		}
		tmp_edges[n1][n2] = 0;
		tmp_edges[n2][n1] = 0;
		tmp_edges[m1][m2] = 0;
		tmp_edges[m2][m1] = 0;

		tmp_edges[n1][m2] = 1;
		tmp_edges[m2][n1] = 1;
		tmp_edges[m1][n2] = 1;
		tmp_edges[n2][m1] = 1;

		temp1 = links + num_links;
		temp2 = links + num_links - 1;
		if ((t1 < num_links - 1 && t2 < num_links - 1) || (t1 >= num_links - 1 && t2 >= num_links - 1)) {
			edge1->first = temp1->first;
			edge1->second = temp1->second;
			edge2->first = temp2->first;
			edge2->second = temp2->second;
		} else {
			if (t1 == num_links) {
				edge2->first = temp2->first;
				edge2->second = temp2->second;
			} else if (t1 == num_links - 1) {
				edge2->first = temp1->first;
				edge2->second = temp1->second;
			}
			if (t2 == num_links) {
				edge1->first = temp2->first;
				edge1->second = temp2->second;
			} else if (t2 == num_links - 1) {
				edge1->first = temp1->first;
				edge1->second = temp1->second;
			}
		}
		temp1->first = n1; temp1->second = m2;
		temp2->first = n2; temp2->second = m1;
	}
}


