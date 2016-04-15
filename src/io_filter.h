/* default parameters for basic machine learning task */
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


extern char *program_invocation_short_name;
/* store program name */
#define ML_PATH "/home/zyp/ML"
//for boolean type
#define FALSE 0
#define TRUE 1
#define False 0
#define True 1
#define F 0
#define T 1

//array macros
#define MAX_NUM_SIZE 20
#define MAX_FILE_NAME_LENGTH 500
#define SINGLE_PRAR_LENGTH 50
#define COMMAND_LINE_STRING_LENGTH 1000

#define MAX_LINE_SIZE 50000
#define MAX_NAME_SIZE 100


#define DELIMITER '\t'

//macro as functions
#define MAX(A, B) A<B?B:A
#define MIN(A, B) A<B?A:B
typedef _Bool boolean;

//for number type
typedef double NUM;
typedef int COUNT;

typedef struct _I2F{// I2F
	COUNT i;
	NUM f;
} I2F;
typedef struct _IJ2F{// I2F
	COUNT i;
	COUNT j;
	NUM f;
} IJ2F;
typedef struct _MN2F{// I2F
	char m[15];
	char n[15];
	NUM f;
} MN2F;

typedef struct _IJ2BPF{// IJ2PF
	COUNT i;
	COUNT j;
	boolean *pf;
} IJ2BPF;
typedef struct _NUM_PF {
	char name[MAX_NAME_SIZE]; 
	NUM *pf;
} NUM_PF;
typedef struct _COUNT_PF {
	char name[MAX_NAME_SIZE]; 
	COUNT *pf;
} COUNT_PF;
typedef struct _BIN_PF {
	char name[MAX_NAME_SIZE]; 
	boolean *pf;
} BIN_PF;
typedef struct _NUM_MAT {
	char **row_name; 
	char **col_name; 
	NUM **mat;
} NUM_MAT;

typedef struct _BIN_ARR{
	int n;
	int *v;
} BIN_ARR;
typedef struct _INT_ARR{
	int nv;
	int nc;
	int *v;
	int *c;
} INT_ARR;

COUNT wc_l(char *file){
	FILE *fp;
	char line[MAX_LINE_SIZE];
	COUNT line_num;

	if ((fp = fopen(file, "r")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
			program_invocation_short_name, file, strerror(errno) );
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
int wc(char *file, COUNT *row_num_ptr, COUNT *col_num_ptr){
	FILE *fp;
	char line[MAX_LINE_SIZE];
	char *token;

	if ((fp = fopen(file, "r")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
			program_invocation_short_name, file, strerror(errno) );
		return -1;
	}
	(*row_num_ptr)=0;
	(*col_num_ptr)=1;
	while (fgets(line, sizeof(line), fp) != NULL 
			&& line[0] != '#' 
			&& line[0] != '^' 
			&& line[0] != '!'){
		if ((*row_num_ptr) == 0) { 
			token = strtok(line ,"\t"); 
			while ( (token = strtok(NULL, "\t\n")) != NULL ) {
				(*col_num_ptr)++;
			}
		}
		(*row_num_ptr)++;
	}
	fclose(fp);
	return 0;
}
void read_str_list (FILE *fp, COUNT row_num, char ***strs) {
	char line[MAX_LINE_SIZE]; 
	char *token; 
	int ri, ri2, col_i;
	*strs = 
		(char **) malloc (row_num * sizeof(char *)); 
	ri = 0; 
	while (fgets(line, sizeof(line), fp) != NULL){
		if ((token = strtok(line, "\t\n")) == NULL)
			break; 
		(*strs)[ri] = 
			( char* )malloc( MAX_NAME_SIZE * sizeof(char) ); 
		strncpy( (*strs)[ri], token, MAX_NAME_SIZE );
		ri++;
	}
}
int read_int_arr_file (char *file, int **arr_ptr, int ln) {
	char line[MAX_LINE_SIZE];
	char *token; 
	int i,j;
	FILE *fp;
	if ((fp = fopen(file, "r")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
			program_invocation_short_name, file, strerror(errno) );
		return -1;
	}
	*arr_ptr = 
		(int *) calloc (ln, sizeof(int)); 
	// skip the first line
	while (fgets(line, sizeof(line), fp) != NULL){
		
		if ((token = strtok(line, "\t")) == NULL)
			break; 
		i = atoi(token);
		if ((token = strtok(NULL, "\t\n")) == NULL)
			break; 
		j = atoi(token);
		(*arr_ptr)[i] = j;
	}
	return 0;
}
int read_num_arrs_indexed (FILE *fp, COUNT mat_rn, COUNT mat_cn, 
	NUM ***arrs_ptr, boolean **index_ptr) {
	char line[MAX_LINE_SIZE]; 
	char *token; 
	int pf_i, pf_i2, col_i;
	*arrs_ptr = 
		(NUM **) malloc (mat_rn * sizeof(NUM *)); 
	*index_ptr = 
		(boolean *) calloc (mat_rn, sizeof(boolean));
	// skip the first line
	while (fgets(line, sizeof(line), fp) != NULL){
		if ((token = strtok(line, "\t")) == NULL)
			break; 
		pf_i = atoi(token);
//		printf("%d\n",pf_i);
		(*index_ptr)[pf_i] = 1;
		(*arrs_ptr)[pf_i] =
			( NUM* )malloc( mat_cn*sizeof(NUM) );
		for(col_i=0; col_i<mat_cn; col_i++){ 
			token = strtok(NULL, "\t\n");
			(*arrs_ptr)[pf_i][col_i] = atof(token);
		}
	}
	return 0;
}
int read_num_pfs (FILE *fp, COUNT row_num, COUNT col_num, 
	NUM_PF **num_pfs_ptr) {
	char line[MAX_LINE_SIZE]; 
	char *token; 
	int pf_i, pf_i2, col_i;
	*num_pfs_ptr = 
		(NUM_PF*) malloc (row_num * sizeof(NUM_PF)); 
	pf_i = 0; 
	while (fgets(line, sizeof(line), fp) != NULL){

		if ((token = strtok(line, "\t")) == NULL)
			break; 
		strncpy((*num_pfs_ptr + pf_i)->name, token, MAX_NAME_SIZE);
		(*num_pfs_ptr + pf_i)->pf = 
			( NUM* )malloc( (col_num-1)*sizeof(NUM) ); 
		for(col_i=0; col_i<col_num-1; col_i++){ 
			token = strtok(NULL, "\t\n");
			(*num_pfs_ptr + pf_i)->pf[col_i] = atof(token);
		}
		pf_i++;
	}
	return 0;
}
int read_bin_pfs (FILE *fp, COUNT row_num, COUNT col_num, 
	BIN_PF **bin_pfs_ptr) {
	char line[MAX_LINE_SIZE]; 
	char *token; 
	int pf_i, pf_i2, col_i;
	*bin_pfs_ptr = 
		(BIN_PF*) malloc (row_num * sizeof(BIN_PF)); 
	// skip the first line
	pf_i = 0; 
	while (fgets(line, sizeof(line), fp) != NULL){

		if ((token = strtok(line, "\t")) == NULL)
			break; 
		strncpy((*bin_pfs_ptr + pf_i)->name, token, MAX_NAME_SIZE);
		(*bin_pfs_ptr + pf_i)->pf = 
			( boolean* )malloc( (col_num-1)*sizeof(boolean) ); 
		for(col_i=0; col_i<col_num-1; col_i++){ 
			token = strtok(NULL, "\t\n");
			if(!atoi(token)){
				(*bin_pfs_ptr + pf_i)->pf[col_i] = TRUE;
			}
			else{
				(*bin_pfs_ptr + pf_i)->pf[col_i] = FALSE;
			}
		}
//		(*num_pfs_ptr + pf_i)->dim = col_num-1;
		pf_i++;
	}
	return 0;
}
void read_num_mat (FILE *fp, COUNT row_num, COUNT col_num, 
	NUM ***mat) {
	char line[MAX_LINE_SIZE]; 
	char *token; 
	int i, j, col_i;
	NUM val;
	*mat = 
		(NUM**) malloc (row_num * sizeof(NUM*)); 
	for(i=0; i<row_num; i++)
		(*mat)[i] = 
			(NUM*) malloc (col_num * sizeof(NUM)); 

	for(i=0; i<row_num; i++)
		for(i=0; i<col_num; i++)
			(*mat)[i][j]=0.0;
	while (fgets(line, sizeof(line), fp) != NULL ){

		if ((token = strtok(line, "\t")) == NULL)
			break; 
		i = atoi(token);
		if ((token = strtok(NULL, "\t")) == NULL)
			break; 
		j = atoi(token);
		if ((token = strtok(NULL, "\t\n")) == NULL)
			break; 
		val = atof(token);
		(*mat)[i][j]=val;
	}
}
void write_count_entry (FILE *fp, COUNT row_num, COUNT col_num, 
	COUNT **mat) {
	COUNT i,j;
	for(i=0; i<row_num; i++){
		for(j=0; j<col_num; j++)
			if(mat[i][j])
				fprintf(fp,"%d\t%d\t%d\n",i,j,mat[i][j]);
	}
}
void write_num_entry (FILE *fp, COUNT row_num, COUNT col_num, 
	NUM **mat) {
	COUNT i,j;
	for(i=0; i<row_num; i++){
		for(j=0; j<col_num; j++)
			if(mat[i][j])
				fprintf(fp,"%d\t%d\t%f\n",i,j,mat[i][j]);
	}
}
void write_num_pfs (FILE *fp, COUNT row_num, COUNT col_num, 
	NUM_PF *num_pfs) {
	COUNT i,j;
	for(i=0; i<row_num; i++){
		fprintf(fp,"%s",num_pfs[i].name);
		for(j=0; j<col_num; j++)
			fprintf(fp,"\t%f",num_pfs[i].pf[j]);
		fprintf(fp,"\n");
	}
}
void write_count_pfs (FILE *fp, COUNT row_num, COUNT col_num, 
	COUNT_PF *cnt_pfs) {
	COUNT i,j;
	for(i=0; i<row_num; i++){
		fprintf(fp,"%d_%s",i,cnt_pfs[i].name);
		for(j=0; j<col_num; j++)
			fprintf(fp,"\t%d",cnt_pfs[i].pf[j]);
		fprintf(fp,"\n");
	}
}





char **store_para(int argc, char *argv[]);
char **split(char *str, char del[2], int *len_ptr);
int compare_num_u (const void *a, const void *b)
{
	return *(NUM*)a > *(NUM*)b ? 1 :-1;
}


int compare_i2f_d (const void *a, const void *b)
{
	return ((I2F *)a)->f < ((I2F *)b)->f ? 1 :-1;
}
void sort_i2f_d(I2F **i2f_arr_ptr, COUNT len){
	qsort((*i2f_arr_ptr), len, sizeof(I2F), compare_i2f_d);
}
void print_i2f_d(I2F **i2f_arr_ptr, COUNT len){
	int i;
	for(i=0;i<len;i++)
		printf("%d\t%f\n",(*i2f_arr_ptr)[i].i,(*i2f_arr_ptr)[i].f);
}



void read_ij2f (FILE *fp, COUNT row_num, 
	IJ2F **ij2fs) {
	char line[MAX_LINE_SIZE]; 
	char *token; 
	int row_i;
	NUM val;
	*ij2fs = 
		(IJ2F*) malloc (row_num * sizeof(IJ2F)); 
	row_i=0;
	while (fgets(line, sizeof(line), fp) != NULL ){

		if ((token = strtok(line, "\t")) == NULL)
			break; 
		(*ij2fs)[row_i].i = atoi(token);
		if ((token = strtok(NULL, "\t\n")) == NULL)
			break; 
		(*ij2fs)[row_i].j = atoi(token);
		if ((token = strtok(NULL, "\t\n")) == NULL){
			(*ij2fs)[row_i].f = 1.0;
			row_i++;
		}
		else{
			(*ij2fs)[row_i].f = atof(token);
			row_i++;
		}
	}
}
int compare_ij2f_d (const void *a, const void *b)
{
	return ((IJ2F *)a)->f < ((IJ2F *)b)->f ? 1 :-1;
}
int compare_ij2f_i (const void *a, const void *b)
{
	return ((IJ2F *)a)->f > ((IJ2F *)b)->f ? 1 :-1;
}
int compare_ij2f_ij (const void *a, const void *b)
{
	if (((IJ2F *)a)->i > ((IJ2F *)b)->i)
		return 1;
	else if (((IJ2F *)a)->i < ((IJ2F *)b)->i)
		return -1;
	else if (((IJ2F *)a)->j > ((IJ2F *)b)->j)
		return 1;
	else
		return -1;
}
char eq_ij2f(IJ2F *a, IJ2F *b) 
{
	if(a->i == b->i 
		&& a->j == b->j)
		return 1;
	return 0;
}
void sort_ij2f_d(IJ2F **ij2f_arr_ptr, COUNT len){
	qsort((*ij2f_arr_ptr), len, sizeof(IJ2F), compare_ij2f_d);
}
void sort_ij2f_i(IJ2F **ij2f_arr_ptr, COUNT len){
	qsort((*ij2f_arr_ptr), len, sizeof(IJ2F), compare_ij2f_i);
}
void sort_ij2f_ij(IJ2F **ij2f_arr_ptr, COUNT len){
	qsort((*ij2f_arr_ptr), len, sizeof(IJ2F), compare_ij2f_ij);
}
void print_ij2f(IJ2F *ij2f_arr, COUNT len){
	int i;
	for(i=0;i<len;i++)
		printf("%d\t%d\t%.3f\n",ij2f_arr[i].i,
			ij2f_arr[i].j,
			ij2f_arr[i].f);
}
void print_file_ij2f(char *file, IJ2F *ij2f_arr, COUNT len){
	int i;
	FILE *fp;

	if ((fp = fopen(file, "w")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
			program_invocation_short_name, file, strerror(errno) );
		return;
	}
	for(i=0;i<len;i++)
		fprintf(fp,"%d\t%d\t%.3f\n",ij2f_arr[i].i,
			ij2f_arr[i].j,
			ij2f_arr[i].f);
	fclose(fp);
}
void print_fp_ij2f_unique(FILE *fp, IJ2F *ij2f_arr, COUNT len){
	int i;
	IJ2F dump;

	dump.i = ij2f_arr[0].i;
	dump.j = ij2f_arr[0].j;
	dump.f = ij2f_arr[0].f;
	for(i=1;i<len;i++){
		if(eq_ij2f(&dump,&ij2f_arr[i])){
			dump.f += ij2f_arr[i].f;
		}
		else{
			fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
			dump.i = ij2f_arr[i].i;
			dump.j = ij2f_arr[i].j;
			dump.f = ij2f_arr[i].f;
		}
	}
	fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
}
void print_file_ij2f_unique(char *file, IJ2F *ij2f_arr, COUNT len){
	int i;
	FILE *fp;
	IJ2F dump;

	if ((fp = fopen(file, "w")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
			program_invocation_short_name, file, strerror(errno) );
		return;
	}
	dump.i = ij2f_arr[0].i;
	dump.j = ij2f_arr[0].j;
	dump.f = ij2f_arr[0].f;
	for(i=1;i<len;i++){
		if(eq_ij2f(&dump,&ij2f_arr[i])){
			dump.f += ij2f_arr[i].f;
		}
		else{
			fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
			dump.i = ij2f_arr[i].i;
			dump.j = ij2f_arr[i].j;
			dump.f = ij2f_arr[i].f;
		}
	}
	fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
	fclose(fp);
}
void print_fp_ij2f_unique_max(FILE *fp, IJ2F *ij2f_arr, COUNT len){
	int i;
	IJ2F dump;

	dump.i = ij2f_arr[0].i;
	dump.j = ij2f_arr[0].j;
	dump.f = ij2f_arr[0].f;
	for(i=1;i<len;i++){
		if(eq_ij2f(&dump,&ij2f_arr[i])){
			if(ij2f_arr[i].f>dump.f)
				dump.f = ij2f_arr[i].f;
		}
		else{
			fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
			dump.i = ij2f_arr[i].i;
			dump.j = ij2f_arr[i].j;
			dump.f = ij2f_arr[i].f;
		}
	}
	fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
}
void print_file_ij2f_unique_max(char *file, IJ2F *ij2f_arr, COUNT len){
	int i;
	FILE *fp;
	IJ2F dump;

	if ((fp = fopen(file, "w")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
			program_invocation_short_name, file, strerror(errno) );
		return;
	}
	dump.i = ij2f_arr[0].i;
	dump.j = ij2f_arr[0].j;
	dump.f = ij2f_arr[0].f;
	for(i=1;i<len;i++){
		if(eq_ij2f(&dump,&ij2f_arr[i])){
			if(ij2f_arr[i].f>dump.f)
				dump.f = ij2f_arr[i].f;
		}
		else{
			fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
			dump.i = ij2f_arr[i].i;
			dump.j = ij2f_arr[i].j;
			dump.f = ij2f_arr[i].f;
		}
	}
	fprintf(fp,"%d\t%d\t%.3f\n",dump.i,dump.j,dump.f);
	fclose(fp);
}


void read_mn2f (FILE *fp, COUNT row_num, 
		MN2F **mn2fs) {
	char line[MAX_LINE_SIZE]; 
	char *token; 
	int row_i;
	NUM val;
	*mn2fs = 
		(MN2F*) malloc (row_num * sizeof(MN2F)); 
	row_i=0;
	while (fgets(line, sizeof(line), fp) != NULL ){

		if ((token = strtok(line, "\t")) == NULL)
			break; 
		strcpy((*mn2fs)[row_i].m, token);
		if ((token = strtok(NULL, "\t")) == NULL)
			break; 
		strcpy((*mn2fs)[row_i].n, token);
		if ((token = strtok(NULL, "\t\n")) == NULL)
			break; 
		(*mn2fs)[row_i].f = atof(token);
		row_i++;
	}
}
int compare_mn2f_d (const void *a, const void *b)
{
	return ((MN2F *)a)->f < ((MN2F *)b)->f ? 1 :-1;
}
int compare_mn2f_i (const void *a, const void *b)
{
	return ((MN2F *)a)->f > ((MN2F *)b)->f ? 1 :-1;
}
void sort_mn2f_d(MN2F **mn2f_arr_ptr, COUNT len){
	qsort((*mn2f_arr_ptr), len, sizeof(MN2F), compare_mn2f_d);
}
void sort_mn2f_i(MN2F **mn2f_arr_ptr, COUNT len){
	qsort((*mn2f_arr_ptr), len, sizeof(MN2F), compare_mn2f_i);
}
void print_mn2f(MN2F *mn2f_arr, COUNT len){
	int i;
	for(i=0;i<len;i++)
		printf("%s\t%s\t%.3f\n",mn2f_arr[i].m,
				mn2f_arr[i].n,
				mn2f_arr[i].f);
}
void print_file_mn2f(char *file, MN2F *mn2f_arr, COUNT len){
	int i;
	FILE *fp;

	if ((fp = fopen(file, "w")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
				program_invocation_short_name, file, strerror(errno) );
		return;
	}
	for(i=0;i<len;i++)
		fprintf(fp,"%s\t%s\t%.3f\n",mn2f_arr[i].m,
				mn2f_arr[i].n,
				mn2f_arr[i].f);
	fclose(fp);
}
void print_file_mn2f_i(char *file, MN2F *mn2f_arr, COUNT len){
	int i;
	FILE *fp;

	if ((fp = fopen(file, "w")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
				program_invocation_short_name, file, strerror(errno) );
		return;
	}
	for(i=0;i<len;i++)
		fprintf(fp,"%s\t%s\t%.0f\n",mn2f_arr[i].m,
				mn2f_arr[i].n,
				mn2f_arr[i].f);
	fclose(fp);
}
void print_file_mn2f_i_cut(char *file, MN2F *mn2f_arr, 
		COUNT len, COUNT cut){
	int i;
	FILE *fp;

	if ((fp = fopen(file, "w")) == NULL){
		fprintf(stderr, "%s: Coudn't open file %s; %s\n",
				program_invocation_short_name, file, strerror(errno) );
		return;
	}
	for(i=0; i<len; i++){
		if( (int) (mn2f_arr[i].f) <= cut)
			break;
		fprintf(fp,"%s\t%s\t%.0f\n",mn2f_arr[i].m,
			mn2f_arr[i].n,
			mn2f_arr[i].f);
	}
	fclose(fp);
}
void ua_init(INT_ARR *ua_ptr){

	ua_ptr->nc=0;
	ua_ptr->nv=0;
	ua_ptr->v=malloc(0);
	ua_ptr->c=malloc(0);
}
void ua_insert(int v, INT_ARR *ua_ptr)
{
	int n,i;
	ua_ptr->nc++;
	n=ua_ptr->nv;
	while (n>0 && v<(ua_ptr->v)[n-1]){
		n--;
	}
	if(n==0 || v>(ua_ptr->v)[n-1]){
		ua_ptr->nv++;
		ua_ptr->v=(int *)realloc(ua_ptr->v, ua_ptr->nv * sizeof(int));
		ua_ptr->c=(int *)realloc(ua_ptr->c, ua_ptr->nv * sizeof(int));
		for(i=ua_ptr->nv-1; i>n; i--){
			(ua_ptr->v)[i] = (ua_ptr->v)[i-1];
			(ua_ptr->c)[i] = (ua_ptr->c)[i-1];
		}
		(ua_ptr->v)[n]=v;
		(ua_ptr->c)[n]=1;
	}
	else{
		(ua_ptr->c)[n-1]++;
	}
}
void ua_insert_range(int stt, int end, INT_ARR *ua_ptr)
{
	int n1,n2,i,flag1,flag2;
	ua_ptr->nc++;
	flag1=0;
	flag2=0;
	n2=ua_ptr->nv;
	while (n2>0 && end<(ua_ptr->v)[n2-1])
		n2--;
	n1=n2;
	while (n1>0 && stt<(ua_ptr->v)[n1-1])
		n1--;
	if (n1!=0 && stt==(ua_ptr->v)[n1-1])
		flag1=1;
	if (n2!=0 && end==(ua_ptr->v)[n2-1])
		flag2=1;
	if(ua_ptr->nv==0 || (!flag1 && !flag2)){
		ua_ptr->nv+=2;
		ua_ptr->v=(int *)realloc(ua_ptr->v, ua_ptr->nv * sizeof(int));
		ua_ptr->c=(int *)realloc(ua_ptr->c, ua_ptr->nv * sizeof(int));

		for(i=ua_ptr->nv-1; i>n2+1; i--){
			(ua_ptr->v)[i] = (ua_ptr->v)[i-2];
			(ua_ptr->c)[i] = (ua_ptr->c)[i-2];
		}
		(ua_ptr->v)[i]=end;
		(ua_ptr->c)[i]=(ua_ptr->c)[i-2];
		for(i=n2; i>n1; i--){
			(ua_ptr->v)[i] = (ua_ptr->v)[i-1];
			(ua_ptr->c)[i] = (ua_ptr->c)[i-1];
			(ua_ptr->c)[i] ++;
		}
		(ua_ptr->v)[i]=stt;
		(ua_ptr->c)[i]=(ua_ptr->c)[i-1];
		(ua_ptr->c)[i]++;
	}
	else if(flag1 && !flag2){
		ua_ptr->nv++;
		ua_ptr->v=(int *)realloc(ua_ptr->v, ua_ptr->nv * sizeof(int));
		ua_ptr->c=(int *)realloc(ua_ptr->c, ua_ptr->nv * sizeof(int));

		for(i=ua_ptr->nv-1; i>n2; i--){
			(ua_ptr->v)[i] = (ua_ptr->v)[i-1];
			(ua_ptr->c)[i] = (ua_ptr->c)[i-1];
		}
		(ua_ptr->v)[i]=end;
		(ua_ptr->c)[i]=(ua_ptr->c)[i-1];
		for(i=n2-1; i>=n1-1; i--){
			(ua_ptr->c)[i] ++;
		}
	}
	else if(flag2 && !flag1){
		ua_ptr->nv++;
		ua_ptr->v=(int *)realloc(ua_ptr->v, ua_ptr->nv * sizeof(int));
		ua_ptr->c=(int *)realloc(ua_ptr->c, ua_ptr->nv * sizeof(int));

		for(i=ua_ptr->nv-1; i>n2-1; i--){
			(ua_ptr->v)[i] = (ua_ptr->v)[i-1];
			(ua_ptr->c)[i] = (ua_ptr->c)[i-1];
		}
		for(i=n2-1; i>n1+1; i--){
			(ua_ptr->c)[i]++;
		}
		printf("%d\t%d\t%d\n",n1,n2,i);
		(ua_ptr->v)[i] = stt;
		(ua_ptr->c)[i] = (ua_ptr->c)[i-1];
		(ua_ptr->c)[i] ++;

	}
	else{
		for(i=n2-2; i>=n1-1; i--){
			(ua_ptr->c)[i]++;
		}
	}
}
void ua_print(INT_ARR ua){
	int i;
	printf("nv: %d\nnc: %d\n",ua.nv,ua.nc);
	printf("v:");
	for(i=0; i<ua.nv; i++){
		printf("\t%d",ua.v[i]);
	}
	printf("\n");
	printf("c:");
	for(i=0; i<ua.nv; i++){
		printf("\t%d",ua.c[i]);
	}
	printf("\n");
}
void ua_destory(INT_ARR *ua_ptr){
	free(ua_ptr->v);
	free(ua_ptr->c);
}

#endif /* _io */
